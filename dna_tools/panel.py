"""
panel.py
========

Kontrola pokrytia CIELOVEHO panelu lokusov (napr. ~1000 pozicii, ktore pouziva
platforma na vypocet genetickych predispozicii) v danom genotypovom subore.

Namiesto imputacie celeho cipu staci overit, ci su tvoje konkretne lokusy
pritomne a zavolane — a cielene imputovat len tie, ktore chybaju.

Vstupny panel (flexibilny parser):
  - obycajny textovy zoznam: jeden rsID na riadok
  - CSV/TSV s hlavickou: rozpozna stlpce rsid / chrom / pos / effect_allele
    (synonyma: rs, snp; chr, chromosome; position, bp; risk_allele, ea, ...)

Match logika (v poradi):
  1. podla rsID (ak panel ma rsID a najde sa v subore)
  2. podla chr:pos (ak panel ma suradnice) — zachytava aj rozdielne rsID verzie
Stav kazdeho lokusu:
  - called   : najdeny, ma platny genotyp
  - no_call  : najdeny, ale genotyp je '--'
  - absent   : v subore vobec nie je (kandidat na imputaciu)
"""

from __future__ import annotations

import csv as _csv
from dataclasses import dataclass, field
from typing import List, Optional, Union

import pandas as pd

from .common import GenotypeData, NO_CALL, normalize_chromosome
from .parsers import parse_file

_RS = {"rsid", "rs", "snp", "snp_id", "marker", "id"}
_CHR = {"chrom", "chr", "chromosome"}
_POS = {"pos", "position", "bp", "coord", "coordinate", "grch37", "b37"}
_EA = {"effect_allele", "risk_allele", "ea", "ra", "allele", "effect"}

import re as _re
_RSID_RE = _re.compile(r"^rs\d+$", _re.I)


def _norm_rsid(s):
    """RS.../Rs... -> rs...  (rsID matchovanie/dohladavanie nema byt case-sensitive).
    i-ID a ine oznacenia ostavaju nezmenene."""
    if s is None:
        return s
    t = str(s).strip()
    if _RSID_RE.match(t):
        return "rs" + t[2:]
    return s


@dataclass
class PanelResult:
    n_panel: int
    called: int
    no_call: int
    absent: int
    matched_by_rsid: int
    matched_by_pos: int
    table: pd.DataFrame            # per-lokus stav
    strand_warnings: int = 0
    enriched: int = 0              # pozicie doplnene z 23andMe referencie (lokalne)
    enriched_online: int = 0       # pozicie doplnene online cez Ensembl

    @property
    def coverage(self) -> float:
        return self.called / self.n_panel if self.n_panel else 0.0

    def summary(self) -> str:
        return (
            "Panel lokusov: %d\n"
            "  zavolane (called) : %d  (%.1f%%)\n"
            "  no-call           : %d\n"
            "  chybajuce (absent): %d   <- kandidati na cielenu imputaciu\n"
            "  match podla rsID  : %d | podla chr:pos: %d\n"
            "  moznE strand/alela nezhody: %d"
            % (self.n_panel, self.called, 100.0 * self.coverage,
               self.no_call, self.absent, self.matched_by_rsid,
               self.matched_by_pos, self.strand_warnings)
        )


def load_panel(path: str) -> pd.DataFrame:
    """Nacita panel do DataFrame so stlpcami: rsid, chromosome, position, effect_allele.

    Chybajuce stlpce budu None. Akceptuje aj obycajny zoznam rsID (jeden na riadok).
    """
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        sample = fh.read(4096)
        fh.seek(0)
        lines = [l.rstrip("\n").rstrip("\r") for l in fh if l.strip()]

    if not lines:
        raise ValueError("Prazdny panel subor.")

    # detekuj oddelovac / ci ma hlavicku
    delim = "\t" if "\t" in sample else ("," if "," in sample else None)

    # obycajny zoznam rsID (bez oddelovaca)
    if delim is None:
        rows = [{"rsid": _norm_rsid(l.strip()), "chromosome": None, "position": None, "effect_allele": None}
                for l in lines if not l.startswith("#")]
        return pd.DataFrame(rows)

    reader = list(_csv.reader(lines, delimiter=delim))
    header = [h.strip().strip('"').lower() for h in reader[0]]
    has_header = any(h in _RS | _CHR | _POS | _EA for h in header)

    def col_idx(names):
        for i, h in enumerate(header):
            if h in names:
                return i
        return None

    if has_header:
        i_rs, i_chr, i_pos, i_ea = col_idx(_RS), col_idx(_CHR), col_idx(_POS), col_idx(_EA)
        body = reader[1:]
    else:
        # bez hlavicky: predpokladaj rsid v 1. stlpci, volitelne chr,pos
        i_rs, i_chr, i_pos, i_ea = 0, (1 if len(header) > 2 else None), (2 if len(header) > 2 else None), None
        body = reader

    rows = []
    for r in body:
        if not r or (r[0].startswith("#")):
            continue
        def get(i):
            return r[i].strip().strip('"') if (i is not None and i < len(r) and r[i].strip()) else None
        rsid = _norm_rsid(get(i_rs))
        chrom = get(i_chr)
        pos = get(i_pos)
        ea = get(i_ea)
        rows.append({
            "rsid": rsid,
            "chromosome": normalize_chromosome(chrom) if chrom else None,
            "position": int(pos) if (pos and pos.isdigit()) else None,
            "effect_allele": ea.upper() if ea else None,
        })
    return pd.DataFrame(rows)


def _lookup_positions(reference_path: str, wanted_rsids) -> dict:
    """Z 23andMe id_position_reference (tsv/tsv.gz) vytiahne {rsid: (chrom, pos)}
    len pre pozadovane rsID (jeden prechod suborom, pamatovo nenarocne)."""
    import gzip
    wanted = set(str(x) for x in wanted_rsids if x)
    if not wanted:
        return {}
    op = gzip.open(reference_path, "rt", encoding="utf-8", errors="replace") \
        if str(reference_path).endswith(".gz") else open(reference_path, "r", encoding="utf-8", errors="replace")
    out = {}
    with op as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h.strip().lower(): i for i, h in enumerate(header)}
        i_id = idx.get("id", idx.get("rsid", 0))
        i_chr = idx.get("chromosome", idx.get("chr", 1))
        i_pos = idx.get("position", idx.get("pos", 2))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_id, i_chr, i_pos):
                continue
            rs = f[i_id].strip()
            if rs in wanted and rs not in out:
                out[rs] = (f[i_chr].strip(), f[i_pos].strip())
                if len(out) == len(wanted):
                    break
    return out


def _lookup_positions_ensembl(rsids, timeout: int = 20) -> dict:
    """Online dohladanie GRCh37 pozicii cez Ensembl REST (davkovo, POST).
    Vrati {rsid: (chrom, pos)}. Pri sietovej chybe vrati co sa podarilo."""
    import json
    import urllib.request
    rs = [str(r) for r in rsids if str(r).startswith("rs")]  # Ensembl pozna len rs...
    if not rs:
        return {}
    canon = set(str(i) for i in range(1, 23)) | {"X", "Y", "MT"}
    url = "https://grch37.rest.ensembl.org/variation/homo_sapiens"
    out = {}
    for i in range(0, len(rs), 190):                    # limit ~200/poziadavka
        body = json.dumps({"ids": rs[i:i + 190]}).encode()
        req = urllib.request.Request(
            url, data=body, method="POST",
            headers={"Content-Type": "application/json", "Accept": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = json.loads(resp.read().decode())
        except Exception:
            continue
        for rid, info in (data or {}).items():
            if not info:
                continue
            for m in info.get("mappings", []):
                chrom = str(m.get("seq_region_name", ""))
                start = m.get("start")
                if chrom in canon and start:
                    out[rid] = (chrom, str(start))
                    break
    return out


def annotate_by_rsid(gd: GenotypeData, id_reference: Optional[str] = None,
                     online: bool = False) -> int:
    """Doplni chromosome/position pre markery, ktore maju platny rsID, ale
    chybajucu/nekanonicku suradnicu (napr. Illumina/GSA rs-only sondy, ktore
    ostali s chromosome 0, lebo v nazve nie je pozicia — len rs cislo).

    Postup: 1) lokalne z 23andMe id_position_reference (zdarma, offline),
            2) ak online=True, zvysne rsID dohlada cez Ensembl GRCh37 REST.
    Namerane/uz-kanonicke suradnice sa NIKDY neprepisu. Vrati pocet doplnenych.
    """
    df = gd.df
    canon = set(str(i) for i in range(1, 23)) | {"X", "Y", "MT"}
    norm_chr = df["chromosome"].astype(str).map(normalize_chromosome)
    rsids = df["rsid"].astype(str)
    need_mask = (~norm_chr.isin(canon)) & rsids.str.match(r"(?i)^rs\d+$")
    wanted = sorted({r.lower() for r in rsids[need_mask]})
    if not wanted:
        return 0

    mapping: dict = {}
    if id_reference:
        mapping.update(_lookup_positions(id_reference, wanted))
    if online:
        missing = [r for r in wanted if r not in mapping]
        if missing:
            mapping.update(_lookup_positions_ensembl(missing))
    if not mapping:
        return 0

    gd.df = df = df.copy()
    # V object priestore doplnime hodnoty, potom stlpce ZJEDNOTIME na jeden typ.
    # (Ak by sme necha­li mix str+int/Int64, neskorsi to_23andme padne na
    #  astype("Int64") s 'cannot safely cast non-equivalent object to int64'.)
    df["chromosome"] = df["chromosome"].astype(object)
    pos = df["position"].astype(object)
    n = 0
    for i in df.index[need_mask]:
        rs = str(df.at[i, "rsid"]).lower()
        hit = mapping.get(rs)
        if not hit or not hit[1]:
            continue
        nc = normalize_chromosome(hit[0])
        if nc not in canon:
            continue
        df.at[i, "chromosome"] = nc
        pos.at[i] = int(float(hit[1]))
        n += 1
    # homogenny nullable Int64 (neplatne/prazdne -> <NA>); to_23andme to zvlada
    df["position"] = pd.to_numeric(pos, errors="coerce").astype("Int64")
    return n


def check_panel(sample: Union[str, GenotypeData], panel: Union[str, pd.DataFrame],
                id_reference: Optional[str] = None, online: bool = False) -> PanelResult:
    """Overi pokrytie panelu lokusov v genotypovom subore/GenotypeData.

    id_reference : ak zadany a panel nema pozicie, doplni chromosome/position
                   z 23andMe referencie (id_position_reference) podla rsID.
    online       : pre rsID nenajdene v referencii dohlada poziciu online (Ensembl).
    """
    data = parse_file(sample) if isinstance(sample, str) else sample
    pdf = load_panel(panel) if isinstance(panel, str) else panel.copy()

    def _fill(mapping) -> int:
        n = 0
        need = pdf["position"].isna() | pdf["chromosome"].isna()
        for i in pdf.index[need]:
            rs = str(pdf.at[i, "rsid"])
            if rs in mapping and mapping[rs][1]:
                pdf.at[i, "chromosome"] = normalize_chromosome(mapping[rs][0])
                try:
                    pdf.at[i, "position"] = int(float(mapping[rs][1]))
                    n += 1
                except (TypeError, ValueError):
                    pass
        return n

    # 1) lokalne z 23andMe referencie (pre rsID-only panely)
    enriched = 0
    if id_reference:
        need = pdf["position"].isna() | pdf["chromosome"].isna()
        if need.any():
            enriched = _fill(_lookup_positions(id_reference, pdf.loc[need, "rsid"].dropna().astype(str)))

    # 2) online pre zvysne rsID (Ensembl GRCh37)
    enriched_online = 0
    if online:
        need = pdf["position"].isna() | pdf["chromosome"].isna()
        if need.any():
            enriched_online = _fill(_lookup_positions_ensembl(pdf.loc[need, "rsid"].dropna().astype(str)))

    df = data.df.copy()
    df["poskey"] = df["chromosome"].astype(str) + ":" + df["position"].astype("Int64").astype(str)
    by_rsid = dict(zip(df["rsid"].astype(str), df["genotype"]))
    by_pos = dict(zip(df["poskey"], df["genotype"]))

    def _cell(v):
        # None / NaN -> None (chrani pred int(NaN) a str(nan))
        return None if v is None or (isinstance(v, float) and pd.isna(v)) else v

    statuses, gts, matched, sw = [], [], [], 0
    for row in pdf.itertuples(index=False):
        gt = None
        how = ""
        rsid = _cell(getattr(row, "rsid", None))
        chrom = _cell(getattr(row, "chromosome", None))
        pos = _cell(getattr(row, "position", None))
        if rsid is not None and str(rsid) in by_rsid:
            gt = by_rsid[str(rsid)]
            how = "rsid"
        elif chrom is not None and pos is not None:
            pk = "%s:%s" % (chrom, int(float(pos)))
            if pk in by_pos:
                gt = by_pos[pk]
                how = "pos"
        if gt is None:
            statuses.append("absent"); gts.append(""); matched.append("")
        elif gt == NO_CALL:
            statuses.append("no_call"); gts.append(gt); matched.append(how)
        else:
            statuses.append("called"); gts.append(gt); matched.append(how)
            # volitelna kontrola: sedi effect_allele s genotypom?
            ea = getattr(row, "effect_allele", None)
            if ea and ea in ("A", "C", "G", "T") and ea not in gt and gt != NO_CALL:
                # effect alela sa v genotype nevyskytuje — moze byt OK (homoz. ref)
                # ale ak ani komplement nesedi na ziadnu alelu, mozny strand problem
                comp = {"A": "T", "T": "A", "C": "G", "G": "C"}[ea]
                if comp in gt and ea not in gt:
                    sw += 1

    out = pdf.copy()
    out["status"] = statuses
    out["genotype"] = gts
    out["matched_by"] = matched

    return PanelResult(
        n_panel=len(out),
        called=statuses.count("called"),
        no_call=statuses.count("no_call"),
        absent=statuses.count("absent"),
        matched_by_rsid=matched.count("rsid"),
        matched_by_pos=matched.count("pos"),
        table=out,
        strand_warnings=sw,
        enriched=enriched,
        enriched_online=enriched_online,
    )
