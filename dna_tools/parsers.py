"""
parsers.py
==========

Parsery, ktore prevedu kazdy podporovany vstupny format do spolocnej internej
reprezentacie (GenotypeData s DataFrame rsid/chromosome/position/genotype).

Hlavny vstupny bod: parse_file(path, fmt=None). Ak fmt nie je zadany, format
sa autodetekuje cez detect.detect_format().

Poznamky ku genotypom:
  - 23andMe, AncestryDNA, MyHeritage, FTDNA hlasia alely na plus vlakne.
  - VCF: rekonstruujeme alely z REF/ALT + GT pola.
  - Illumina Final Report: alely su typicky na TOP/BOT vlakne; bez manifestu
    (RefStrand) parser prida varovanie. Preferuje stlpce "Allele - Plus".
"""

from __future__ import annotations

import csv
import re
from typing import Optional

import pandas as pd

from .common import (
    GenotypeData,
    NO_CALL,
    CANONICAL_CHROMS,
    normalize_chromosome,
    normalize_genotype,
)
from .detect import FormatType, detect_format, _open_text


# Vzory na cistenie Illumina/GSA nazvov sond.
#   - rsID moze byt "obaleny" prefixom/suffixom: BOT-rs1135675, seq-rs3209663,
#     GSA-rs..., rs111647200_ilmndup1 ...  -> vytiahni holy rs<cislo>.
#   - pozicia moze byt zapisana priamo v nazve: 1:110228436_CNV_GSTM1,
#     22:24301858_CNV_GSTT2B_Ilmndup1, chr7:6026426 ... -> vytiahni chrom+poziciu.
_RS_RE = re.compile(r"(rs\d+)", re.IGNORECASE)
_CHRPOS_RE = re.compile(r"^(?:chr)?(\d{1,2}|X|Y|XY|MT|M)[:_-](\d+)(?![0-9])", re.IGNORECASE)
# Nekotvene vzory pre chr:pos v STREDE nazvu (napr. 'IlmnSeq_10:86912287_IlmnFwd',
# 'IDS-chrX-148584849'). Pouzivaju sa LEN ked stlpce Chr/Position su 0/prazdne.
_CHRPOS_MID_RE = re.compile(r"(?:chr)?(\d{1,2}|X|Y|MT):(\d+)", re.IGNORECASE)
_CHRPOS_MID_DASH_RE = re.compile(r"chr(\d{1,2}|X|Y|MT)-(\d+)", re.IGNORECASE)
_CANON = set(CANONICAL_CHROMS)


def clean_illumina_name(name: str, chrom: str, pos: str):
    """Z Illumina/GSA nazvu sondy odvodi cistejsi (rsid, chromosome, position).

    - Ak nazov obsahuje dbSNP rsID (aj obaleny textom), pouzije sa holy 'rs<cislo>'.
    - Inak, ak je v nazve citatelna 'chrom:pozicia', nazov sa skrati na 'chrom:pos'.
    - Chromozom/poziciu doplni Z NAZVU iba vtedy, ked chybaju v stlpcoch Chr/Position
      (Illumina davá 0 pri nemapovanych CNV/indel sondach) — namerane suradnice
      zo stlpcov maju vzdy prednost.
    Vrati (rsid, chromosome, position).
    """
    name = (name or "").strip()
    col_chrom = normalize_chromosome(chrom)
    col_pos = (pos or "").strip()
    have_chrom = col_chrom in _CANON
    have_pos = col_pos not in ("", "0")

    # 1) suradnice citatelne z nazvu (pouzijeme len na doplnenie chybajucich)
    name_chrom = name_pos = None
    m = _CHRPOS_RE.match(name)
    if m:
        nc = normalize_chromosome(m.group(1))
        if nc in _CANON:
            name_chrom, name_pos = nc, m.group(2)

    # 1b) ak stlpce chybaju A chr:pos nie je na zaciatku nazvu, skus najst chr:pos
    #     kdekolvek v strede nazvu (napr. 'IlmnSeq_10:86912287_IlmnFwd',
    #     'IDS-chrX-148584849'). Namerane stlpce maju vzdy prednost, preto to
    #     robime len ked chyba chrom alebo pozicia a start-anchor nezabral.
    if name_chrom is None and (not have_chrom or not have_pos):
        mm = _CHRPOS_MID_RE.search(name) or _CHRPOS_MID_DASH_RE.search(name)
        if mm:
            nc = normalize_chromosome(mm.group(1))
            if nc in _CANON:
                name_chrom, name_pos = nc, mm.group(2)

    out_chrom = col_chrom if have_chrom else (name_chrom or col_chrom)
    out_pos = col_pos if have_pos else (name_pos or col_pos)

    # 2) rsid — skryte rs cislo ma prednost, inak cista chrom:pos, inak povodny nazov
    rs = _RS_RE.search(name)
    if rs:
        rsid = rs.group(1).lower()
    elif name_chrom and name_pos:
        rsid = "%s:%s" % (name_chrom, name_pos)
    else:
        rsid = name

    return rsid, out_chrom, out_pos


def _iter_noncomment(path: str, comment_prefixes=("#",)):
    """Generator riadkov, ktory preskoci komentarove a prazdne riadky."""
    with _open_text(path) as fh:
        for line in fh:
            s = line.rstrip("\n").rstrip("\r")
            if not s.strip():
                continue
            if any(s.startswith(p) for p in comment_prefixes):
                continue
            yield s


# --- DTC textove formaty ---------------------------------------------------

def parse_23andme(path: str) -> GenotypeData:
    """23andMe: tab-oddelene rsid, chromosome, position, genotype."""
    rows = []
    for s in _iter_noncomment(path):
        parts = s.split("\t")
        if len(parts) < 4:
            continue
        rsid, chrom, pos, gt = parts[0], parts[1], parts[2], parts[3]
        if rsid.lower() == "rsid":
            continue
        rows.append((rsid, normalize_chromosome(chrom), pos, normalize_genotype(gt)))
    df = pd.DataFrame(rows, columns=["rsid", "chromosome", "position", "genotype"])
    return GenotypeData(df, source_format=FormatType.TWENTYTHREEANDME.value)


def parse_ancestrydna(path: str) -> GenotypeData:
    """AncestryDNA: tab-oddelene rsid, chromosome, position, allele1, allele2."""
    rows = []
    for s in _iter_noncomment(path):
        parts = s.split("\t")
        if len(parts) < 5:
            continue
        rsid, chrom, pos, a1, a2 = parts[0], parts[1], parts[2], parts[3], parts[4]
        if rsid.lower() == "rsid":
            continue
        rows.append((rsid, normalize_chromosome(chrom), pos, normalize_genotype(a1, a2)))
    df = pd.DataFrame(rows, columns=["rsid", "chromosome", "position", "genotype"])
    gd = GenotypeData(df, source_format=FormatType.ANCESTRYDNA.value)
    gd.warnings.append("AncestryDNA koduje no-call ako '0' — mapovane na '--'.")
    return gd


def _parse_myheritage_like(path: str, fmt: str) -> GenotypeData:
    """MyHeritage / FTDNA: comma-oddelene, quoted: RSID, CHROMOSOME, POSITION, RESULT."""
    rows = []
    with _open_text(path) as fh:
        reader = csv.reader(fh)
        for parts in reader:
            if not parts or parts[0].startswith("#"):
                continue
            if len(parts) < 4:
                continue
            rsid = parts[0].strip().strip('"')
            if rsid.lower() == "rsid":
                continue
            chrom = parts[1].strip().strip('"')
            pos = parts[2].strip().strip('"')
            res = parts[3].strip().strip('"')
            rows.append((rsid, normalize_chromosome(chrom), pos, normalize_genotype(res)))
    df = pd.DataFrame(rows, columns=["rsid", "chromosome", "position", "genotype"])
    return GenotypeData(df, source_format=fmt)


def parse_myheritage(path: str) -> GenotypeData:
    return _parse_myheritage_like(path, FormatType.MYHERITAGE.value)


def parse_ftdna(path: str) -> GenotypeData:
    return _parse_myheritage_like(path, FormatType.FTDNA.value)


# --- VCF -------------------------------------------------------------------

def parse_vcf(path: str) -> GenotypeData:
    """Parsuje bialelicky VCF (prva vzorka). Haploidne GT (napr. '1') -> 1 znak."""
    rows = []
    n_skipped = 0
    sample_name = None
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                if len(cols) > 9:
                    sample_name = cols[9]
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            alleles = [ref] + alt.split(",")
            if any(len(a) != 1 or a in (".", "*") for a in alleles):
                n_skipped += 1
                continue
            gt = None
            if len(parts) >= 10:
                fmt_keys = parts[8].split(":")
                sample_vals = parts[9].split(":")
                if "GT" in fmt_keys:
                    gt = sample_vals[fmt_keys.index("GT")]
            if gt is None:
                n_skipped += 1
                continue
            gt_norm = gt.replace("|", "/")
            if gt_norm in (".", "./."):
                genotype = NO_CALL
            else:
                idx = gt_norm.split("/")
                try:
                    a1 = alleles[int(idx[0])]
                    if len(idx) > 1:
                        a2 = alleles[int(idx[1])]
                        genotype = normalize_genotype(a1, a2)
                    else:
                        genotype = normalize_genotype(a1)   # haploid -> 1 znak
                except (ValueError, IndexError):
                    genotype = NO_CALL
            rsid = vid if vid and vid != "." else "%s:%s" % (normalize_chromosome(chrom), pos)
            rows.append((rsid, normalize_chromosome(chrom), pos, genotype))
    df = pd.DataFrame(rows, columns=["rsid", "chromosome", "position", "genotype"])
    gd = GenotypeData(df, source_format=FormatType.VCF.value, sample_id=sample_name)
    if n_skipped:
        gd.warnings.append("VCF: preskocenych %d ne-SNP / multialelickych / bez-GT zaznamov." % n_skipped)
    gd.warnings.append("VCF: over, ze suradnice su v cielovom builde (23andMe = GRCh37).")
    return gd


# --- PLINK -----------------------------------------------------------------

def parse_plink(ped_path: str, map_path: Optional[str] = None) -> GenotypeData:
    """Parsuje PLINK .ped/.map (prva vzorka). map_path sa odvodi zamenou pripony."""
    if map_path is None:
        if ped_path.endswith(".ped"):
            map_path = ped_path[:-4] + ".map"
        else:
            raise ValueError("Pre PLINK treba zadat .map subor (map_path).")
    variants = []
    for s in _iter_noncomment(map_path):
        parts = s.split()
        if len(parts) < 4:
            continue
        variants.append((normalize_chromosome(parts[0]), parts[1], parts[3]))
    sample_id = None
    alleles = []
    for s in _iter_noncomment(ped_path):
        parts = s.split()
        if len(parts) < 6:
            continue
        sample_id = parts[1]
        alleles = parts[6:]
        break
    rows = []
    for i, (chrom, rsid, pos) in enumerate(variants):
        a1 = alleles[2 * i] if 2 * i < len(alleles) else "0"
        a2 = alleles[2 * i + 1] if 2 * i + 1 < len(alleles) else "0"
        rows.append((rsid, chrom, pos, normalize_genotype(a1, a2)))
    df = pd.DataFrame(rows, columns=["rsid", "chromosome", "position", "genotype"])
    return GenotypeData(df, source_format=FormatType.PLINK_PED.value, sample_id=sample_id)


# --- Illumina Final Report -------------------------------------------------

def parse_illumina_final_report(path: str) -> GenotypeData:
    """Parsuje Illumina Final Report; preferuje Plus vlakno, inak Top s varovanim."""
    with _open_text(path) as fh:
        lines = [l.rstrip("\n").rstrip("\r") for l in fh]
    data_idx = None
    for i, l in enumerate(lines):
        if l.strip().lower() == "[data]":
            data_idx = i
            break
    header_idx = data_idx + 1 if data_idx is not None else 0
    while header_idx < len(lines) and not lines[header_idx].strip():
        header_idx += 1
    # Illumina Final Report byva tab- ALEBO comma-oddeleny -> rozpoznaj oddelovac
    header_line = lines[header_idx]
    delim = "\t" if "\t" in header_line else ("," if "," in header_line else "\t")
    header = header_line.split(delim)
    col = {name.strip().lower(): j for j, name in enumerate(header)}

    def find(*names):
        for n in names:
            if n in col:
                return col[n]
        return None

    i_snp = find("snp name", "name")
    i_chr = find("chr", "chromosome")
    i_pos = find("position", "mapinfo", "pos")
    plus1, plus2 = find("allele1 - plus"), find("allele2 - plus")
    fwd1, fwd2 = find("allele1 - forward"), find("allele2 - forward")
    top1, top2 = find("allele1 - top"), find("allele2 - top")
    if plus1 is not None and plus2 is not None:
        i_a1, i_a2, used = plus1, plus2, "Plus"
    elif fwd1 is not None and fwd2 is not None:
        i_a1, i_a2, used = fwd1, fwd2, "Forward"
    elif top1 is not None and top2 is not None:
        i_a1, i_a2, used = top1, top2, "Top"
    else:
        raise ValueError("Illumina Final Report: nenasli sa stlpce s alelami.")
    rows = []
    n_rs_extracted = 0     # rsID vytiahnuty z obaleneho nazvu (BOT-rs... -> rs...)
    n_coords_recovered = 0  # chrom/pos doplnene z nazvu (stlpce boli 0/prazdne)
    n_named_by_pos = 0     # marker oznaceny suradnicou (chybal stlpec SNP Name)
    n_unnamed = 0          # ani nazov, ani pouzitelna suradnica
    for l in lines[header_idx + 1:]:
        if not l.strip():
            continue
        parts = l.split(delim)
        if max(i_snp or 0, i_chr or 0, i_pos or 0, i_a1, i_a2) >= len(parts):
            continue
        raw_name = parts[i_snp] if i_snp is not None else ""
        raw_chrom = parts[i_chr] if i_chr is not None else "0"
        raw_pos = parts[i_pos] if i_pos is not None else "0"
        rsid, chrom, pos = clean_illumina_name(raw_name, raw_chrom, raw_pos)
        if rsid in ("", "."):
            # Export bez stlpca 'SNP Name' (napr. Eurofins Orion_v2). Marker
            # oznacime suradnicou — relabel_from_reference mu potom podla
            # (chrom, pozicia) doplni skutocne rsID z 23andMe referencie.
            if chrom in _CANON and pos not in ("", "0"):
                rsid = "%s:%s" % (chrom, pos)
                n_named_by_pos += 1
            else:
                rsid = "."
                n_unnamed += 1
        if rsid != raw_name.strip():
            if rsid.startswith("rs"):
                n_rs_extracted += 1
        if normalize_chromosome(raw_chrom) not in _CANON and chrom in _CANON:
            n_coords_recovered += 1
        gt = normalize_genotype(parts[i_a1], parts[i_a2])
        rows.append((rsid, chrom, pos, gt))
    df = pd.DataFrame(rows, columns=["rsid", "chromosome", "position", "genotype"])
    gd = GenotypeData(df, source_format=FormatType.ILLUMINA_FINAL_REPORT.value)
    if n_rs_extracted:
        gd.warnings.append(
            "Illumina: %d rsID vytiahnutych z obalenych nazvov sond (napr. 'BOT-rs123' -> 'rs123')." % n_rs_extracted)
    if n_coords_recovered:
        gd.warnings.append(
            "Illumina: %d sondam doplnene chrom/pozicia z nazvu (stlpce Chr/Position boli 0)." % n_coords_recovered)
    if i_snp is None:
        gd.warnings.append(
            "Illumina: export nema stlpec 'SNP Name' — %d markerov oznacenych ako 'chrom:pozicia'. "
            "Skutocne rsID doplni az relabel podla 23andMe referencie (relabel_from_reference)."
            % n_named_by_pos)
    if n_unnamed:
        gd.warnings.append(
            "Illumina: %d markerov ostalo bez identifikatora (chyba nazov aj pouzitelna pozicia)." % n_unnamed)
    if used in ("Top", "Forward"):
        gd.warnings.append(
            "Illumina: pouzite '%s' vlakno. 23andMe je na Plus (forward) vlakne. "
            "Treba strand-flip podla RefStrand z manifestu (.bpm) alebo `bcftools +fixref`." % used)
    else:
        gd.warnings.append("Illumina: pouzite Plus vlakno — kompatibilne s 23andMe.")
    return gd


# --- Preznacenie podla 23andMe referencie (podla pozicie) ------------------

def relabel_from_reference(gd: GenotypeData, reference_path: str) -> int:
    """Pre kazdy marker najde jeho (chrom,pos) v 23andMe referencii a prevezme
    tamojsie oznacenie (id) do stlpca rsid. Kde sa pozicia nenajde, rsid ostava
    nezmeneny. Nic sa nevytahuje z povodnych nazvov. Vrati pocet preznacenych."""
    df = gd.df
    chroms = df["chromosome"].astype(str).map(normalize_chromosome)
    positions = pd.to_numeric(df["position"], errors="coerce").astype("Int64")
    want = set()
    for c, p in zip(chroms, positions):
        if p is not pd.NA:
            want.add((c, int(p)))
    if not want:
        return 0
    ref = {}
    with _open_text(reference_path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        idx = {h.strip().lower(): j for j, h in enumerate(head)}
        i_id = idx.get("id", idx.get("rsid", 0))
        i_chr = idx.get("chromosome", idx.get("chr", 1))
        i_pos = idx.get("position", idx.get("pos", 2))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(i_id, i_chr, i_pos):
                continue
            try:
                key = (normalize_chromosome(f[i_chr]), int(f[i_pos]))
            except ValueError:
                continue
            if key in want and key not in ref:
                ref[key] = f[i_id].strip()
    n = 0
    new_rsid = []
    for c, p, old in zip(chroms, positions, df["rsid"]):
        k = (c, int(p)) if p is not pd.NA else None
        if k in ref and ref[k]:
            new_rsid.append(ref[k])
            n += 1
        else:
            new_rsid.append(old)
    gd.df = gd.df.copy()
    gd.df["rsid"] = new_rsid
    return n


# --- Dispatcher ------------------------------------------------------------

_PARSERS = {
    FormatType.TWENTYTHREEANDME: parse_23andme,
    FormatType.ANCESTRYDNA: parse_ancestrydna,
    FormatType.MYHERITAGE: parse_myheritage,
    FormatType.FTDNA: parse_ftdna,
    FormatType.VCF: parse_vcf,
    FormatType.ILLUMINA_FINAL_REPORT: parse_illumina_final_report,
}


def parse_file(path: str, fmt: Optional[FormatType] = None, map_path: Optional[str] = None) -> GenotypeData:
    """Parsuje subor do GenotypeData. Ak fmt nie je zadany, autodetekuje sa."""
    if fmt is None:
        fmt = detect_format(path)
    if fmt in (FormatType.PLINK_PED, FormatType.PLINK_MAP):
        ped = path if path.endswith(".ped") else (path[:-4] + ".ped")
        return parse_plink(ped, map_path)
    parser = _PARSERS.get(fmt)
    if parser is None:
        raise ValueError("Nepodporovany alebo nerozpoznany format: %s" % fmt)
    return parser(path)
