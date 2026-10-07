"""
targeted_impute.py
==================

Cielena imputacia len konkretnych lokusov (napr. panel ~1000 pozicii), nie
celeho cipu. Namiesto imputacie celeho genomu sa okolo cielov vyrezu uzke
okna (default +-250 kb), imputuju sa len tie a z vysledku sa vytiahnu cielove
genotypy so skore spolahlivosti (Beagle DR2).

Tento modul obsahuje CISTO-PYTHONOVE, offline-testovatelne casti:
  - build_windows()            : cielove pozicie -> zluceene regiony (BED)
  - extract_targets_from_vcf() : z imputovaneho VCF vytiahne GT+DR2 pre ciele
  - merge_panel()              : skombinuje 'called' zo vzorky + imputovane ciele

Samotny Beagle beh (ktory potrebuje 1000G referenciu) je v targeted_impute.sh.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pandas as pd

from .common import NO_CALL, normalize_chromosome, normalize_genotype, sorted_genotype, chrom_sort_key


# ---------------------------------------------------------------------------
# 1) Okna okolo cielov
# ---------------------------------------------------------------------------

@dataclass
class Region:
    chrom: str
    start: int
    end: int
    targets: List[int] = field(default_factory=list)   # pozicie cielov v regione

    def __str__(self):
        return f"{self.chrom}:{self.start}-{self.end} ({len(self.targets)} cielov)"


def build_windows(targets: List[Tuple[str, int]], flank: int = 250_000) -> List[Region]:
    """Z cielovych (chrom,pos) vytvori zlucene regiony +-flank.

    Prekryvajuce sa / dotykajuce sa okna na tom istom chromozome sa zlucia,
    aby sa Beagle spustil na co najmensom pocte regionov.
    """
    # normalizuj a zorad
    pts = sorted(
        ((normalize_chromosome(c), int(p)) for c, p in targets),
        key=lambda cp: (chrom_sort_key(cp[0]), cp[1]),
    )
    regions: List[Region] = []
    for chrom, pos in pts:
        s, e = max(1, pos - flank), pos + flank
        if regions and regions[-1].chrom == chrom and s <= regions[-1].end:
            regions[-1].end = max(regions[-1].end, e)
            regions[-1].targets.append(pos)
        else:
            regions.append(Region(chrom, s, e, [pos]))
    return regions


def write_regions_bed(regions: List[Region], path: str) -> str:
    """Zapise regiony do BED (0-based start), pouzitelne pre bcftools/tabix."""
    with open(path, "w") as fh:
        for r in regions:
            fh.write(f"{r.chrom}\t{max(0, r.start-1)}\t{r.end}\n")
    return path


# ---------------------------------------------------------------------------
# 2) Extrakcia cielov z imputovaneho VCF (GT + DR2)
# ---------------------------------------------------------------------------

def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _alleles_to_23andme(a1: str, a2: str, all_alleles) -> str:
    """Genotyp z VCF alel -> 23andMe zapis.
    SNP -> bazy (abecedne, napr. 'CT'). Indel -> I/D konvencia:
    dlhsia alela = insertion 'I', kratsia = deletion 'D' (napr. TC/T -> 'DI')."""
    if any(len(a) != 1 for a in all_alleles):
        mx = max(len(a) for a in all_alleles)
        mn = min(len(a) for a in all_alleles)

        def code(a):
            if mx == mn:
                return "I"
            return "I" if len(a) == mx else "D"

        return "".join(sorted([code(a1), code(a2)]))
    return sorted_genotype(normalize_genotype(a1, a2))


def extract_targets_from_vcf(vcf_path: str,
                             targets: Dict[Tuple[str, int], str],
                             min_gp: float = 0.9) -> pd.DataFrame:
    """Z imputovaneho VCF vytiahne genotypy na cielovych poziciach.

    targets : {(chrom,pos): rsid} — ktore pozicie a pod akym rsID ich chceme.
    min_gp  : prah spolahlivosti. Pri imputacii JEDNEJ vzorky je spravna miera
              istoty max(GP) (najvyssia pravdepodobnost genotypu), NIE DR2
              (ktore je pri 1 vzorke degenerovane). Genotypy s max_gp < min_gp
              sa oznacia priznakom low_conf.
    Vrati DataFrame: rsid, chromosome, position, genotype, max_gp, dr2, low_conf.
    """
    # index podla pozicie aj podla rsID (rsID zachyti indely s 1bp posunom kotvenia)
    want = {(normalize_chromosome(c), int(p)): rs for (c, p), rs in targets.items()}
    want_rsid = {str(rs).lower(): (normalize_chromosome(c), int(p))
                 for (c, p), rs in targets.items() if str(rs).lower().startswith("rs")}
    matched = set()   # ktore cielove (chrom,pos) uz mame
    rows = []
    with _open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos = normalize_chromosome(f[0]), int(f[1])
            vcf_rsid = f[2].lower()
            # 1) match podla pozicie, inak 2) podla rsID (napr. indel o 1bp vedla)
            if (chrom, pos) in want:
                tkey, trsid = (chrom, pos), want[(chrom, pos)]
            elif vcf_rsid in want_rsid and want_rsid[vcf_rsid] not in matched:
                tkey = want_rsid[vcf_rsid]
                trsid = want.get(tkey, f[2])
            else:
                continue
            if tkey in matched:
                continue
            matched.add(tkey)
            ref, alt = f[3], f[4].split(",")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            try:
                dr2 = float(str(info.get("DR2", info.get("R2", "nan"))).replace(",", "."))
            except ValueError:
                dr2 = float("nan")
            alleles = [ref] + alt
            fmt_keys = f[8].split(":")
            sample_vals = f[9].split(":")
            fv = dict(zip(fmt_keys, sample_vals))
            # max(GP): najvyssia pravdepodobnost genotypu (spolahlivost pre 1 vzorku)
            max_gp = float("nan")
            if "GP" in fv:
                try:
                    gps = [float(x.replace(",", ".")) for x in fv["GP"].split(",")]
                    if gps:
                        max_gp = max(gps)
                except ValueError:
                    pass
            gt = fv.get("GT", ".").replace("|", "/")
            if gt in (".", "./."):
                genotype = NO_CALL
            else:
                idx = gt.split("/")
                try:
                    a1 = alleles[int(idx[0])]
                    a2 = alleles[int(idx[1])] if len(idx) > 1 else a1
                    genotype = _alleles_to_23andme(a1, a2, alleles)
                except (ValueError, IndexError):
                    genotype = NO_CALL
            low_conf = (max_gp == max_gp and max_gp < min_gp)
            rows.append({
                "rsid": trsid, "chromosome": tkey[0], "position": tkey[1],
                "genotype": genotype, "max_gp": max_gp, "dr2": dr2,
                "low_conf": low_conf,
            })
    return pd.DataFrame(rows, columns=["rsid", "chromosome", "position",
                                       "genotype", "max_gp", "dr2", "low_conf"])


# ---------------------------------------------------------------------------
# 3) Zlucenie: 'called' zo vzorky + imputovane ciele -> finalny panel
# ---------------------------------------------------------------------------

def merge_panel(panel_status_csv: str,
                imputed: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Skombinuje vysledok panel-check (called/absent) s imputovanymi genotypmi.

    Vrati finalny panel: rsid, chromosome, position, genotype, source, confidence, dr2.
    source     = 'measured' | 'imputed' | 'missing'
    confidence = max(GP) pre imputovane (1.0 = ista), prazdne pre measured/missing.
    """
    base = pd.read_csv(panel_status_csv, dtype=str)
    imp_map = {}
    if imputed is not None and len(imputed):
        for r in imputed.itertuples(index=False):
            imp_map[str(r.rsid)] = (r.genotype, getattr(r, "max_gp", ""), getattr(r, "dr2", ""))

    out = []
    for r in base.itertuples(index=False):
        rsid = str(r.rsid)
        status = getattr(r, "status", "")
        if status == "called":
            out.append((rsid, r.chromosome, r.position, r.genotype, "measured", "", ""))
        elif rsid in imp_map:
            gt, gp, dr2 = imp_map[rsid]
            conf = round(gp, 3) if isinstance(gp, float) and gp == gp else ""
            out.append((rsid, r.chromosome, r.position, gt, "imputed", conf, dr2))
        else:
            out.append((rsid, r.chromosome, r.position, NO_CALL, "missing", "", ""))
    return pd.DataFrame(out, columns=["rsid", "chromosome", "position",
                                      "genotype", "source", "confidence", "dr2"])
