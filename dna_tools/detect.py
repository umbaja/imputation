"""
detect.py
=========

Automaticka detekcia formatu surovych genetickych dat.

DOLEZITE: detekcia je STRUKTURNA (podla oddelovaca a stlpcov dat), nie podla
volnych klucovych slov v komentaroch. Preto sa napr. prevedeny 23andMe subor,
ktory v hlavicke spomina zdrojovy format, spravne rozpozna ako 23andMe (tab,
4 stlpce), a nie ako MyHeritage.

Priorita:
  1. autoritativne hlavicky: ##fileformat=VCF / ##fileformat=MyHeritage
  2. Illumina [Header]/[Data] sekcie
  3. struktura datovej casti (tab vs comma, pocet stlpcov, allele1/allele2)
  4. komentarove klucove slova (len ako posledna moznost)
"""

from __future__ import annotations

import gzip
import os
from enum import Enum
from typing import List, Tuple


class FormatType(str, Enum):
    TWENTYTHREEANDME = "23andme"
    ANCESTRYDNA = "ancestrydna"
    MYHERITAGE = "myheritage"
    FTDNA = "ftdna"
    VCF = "vcf"
    PLINK_PED = "plink_ped"
    PLINK_MAP = "plink_map"
    ILLUMINA_FINAL_REPORT = "illumina_final_report"
    UNKNOWN = "unknown"


def _open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _read_head(path: str, n_lines: int = 80) -> List[str]:
    lines: List[str] = []
    with _open_text(path) as fh:
        for line in fh:
            lines.append(line.rstrip("\n").rstrip("\r"))
            if len(lines) >= n_lines:
                break
    return lines


def _first_data_line(head: List[str]):
    """Vrati (header_or_data_line, is_header) — prvy neprazdny riadok, ktory nie je cisty komentar."""
    for l in head:
        if not l.strip():
            continue
        if l.startswith("##"):
            continue
        if l.startswith("#"):
            # 23andMe ma stlpcovu hlavicku ako komentar: "# rsid chromosome position genotype"
            low = l.lower()
            if "rsid" in low and ("genotype" in low or "position" in low):
                return l.lstrip("#").strip(), True
            continue
        return l, False
    return None, False


def _illumina_data_header(head: List[str]) -> List[str]:
    """Vrati stlpcovu hlavicku Illumina Final Reportu (prvy neprazdny riadok
    za [Data]) ako mala pismena. Sekcia [Header] a sekcia [Data] mozu mat RUZNE
    oddelovace (Eurofins: [Header] tabmi, [Data] ciarkami), preto oddelovac
    urcujeme az z toho konkretneho riadku."""
    for i, l in enumerate(head):
        if l.strip().lower() != "[data]":
            continue
        for l2 in head[i + 1:]:
            if not l2.strip():
                continue
            delim = "\t" if "\t" in l2 else ","
            return [c.strip().strip('"').lower() for c in l2.split(delim)]
        return []
    return []


def detect_format_verbose(path: str) -> Tuple[FormatType, str]:
    ext = os.path.splitext(path.replace(".gz", ""))[1].lower()
    if ext == ".ped":
        return FormatType.PLINK_PED, "pripona .ped (PLINK genotypy)"
    if ext == ".map":
        return FormatType.PLINK_MAP, "pripona .map (PLINK pozicie)"

    head = _read_head(path)
    if not head:
        return FormatType.UNKNOWN, "prazdny subor"
    lower = "\n".join(head).lower()

    # --- 1. autoritativne hlavicky ---
    if head[0].startswith("##fileformat=VCF") or "##fileformat=vcf" in lower:
        return FormatType.VCF, "hlavicka ##fileformat=VCF"
    if any(l.startswith("#CHROM\tPOS\tID\tREF\tALT") for l in head):
        return FormatType.VCF, "VCF stlpcova hlavicka #CHROM POS ID REF ALT"
    if any(l.lower().startswith("##fileformat=myheritage") for l in head):
        return FormatType.MYHERITAGE, "autoritativna hlavicka ##fileformat=MyHeritage"

    # --- 2. Illumina ---
    # Rozpoznavame STRUKTURALNE: podla stlpcovej hlavicky hned za [Data], nie
    # podla pritomnosti konkretneho nazvu stlpca. Rozne linky (GSGT, Eurofins
    # Orion_v2, ...) exportuju rozne sady stlpcov — jedine spolocne je dvojica
    # Allele1/Allele2 - <vlakno>. Stlpec 'SNP Name' moze aj chybat.
    if "[data]" in lower or "[header]" in lower:
        hdr_cols = _illumina_data_header(head)
        has_a1 = any(c.startswith("allele1 - ") for c in hdr_cols)
        has_a2 = any(c.startswith("allele2 - ") for c in hdr_cols)
        if has_a1 and has_a2:
            alle = ", ".join(c for c in hdr_cols if c.startswith("allele"))
            extra = "" if "snp name" in hdr_cols else " (bez stlpca 'SNP Name')"
            return (FormatType.ILLUMINA_FINAL_REPORT,
                    "sekcie [Header]/[Data] a stlpce %s%s" % (alle, extra))
        if "snp name" in hdr_cols or "gc score" in hdr_cols:
            return (FormatType.ILLUMINA_FINAL_REPORT,
                    "sekcie [Header]/[Data] a stlpec 'SNP Name' / 'GC Score'")
        if "allele1 - top" in lower or "snp name" in lower or "gc score" in lower:
            return (FormatType.ILLUMINA_FINAL_REPORT,
                    "sekcie [Header]/[Data] a stlpce 'SNP Name' / 'Allele1 - Top'")

    # --- 3. struktura datovej casti ---
    line, is_header = _first_data_line(head)
    if line is not None:
        low = line.lower()
        tab = line.split("\t")
        com = [c.strip().strip('"') for c in line.split(",")]

        # comma-oddelene -> MyHeritage / FTDNA
        if len(com) >= 4 and (("rsid" in low and "result" in low) or com[0].startswith("rs")):
            if "familytreedna" in lower or "ftdna" in lower:
                return FormatType.FTDNA, "comma-oddelene data (FamilyTreeDNA podla komentara)"
            return FormatType.MYHERITAGE, "comma-oddelene data RSID,CHROMOSOME,POSITION,RESULT (MyHeritage/FTDNA)"

        # tab-oddelene
        if len(tab) >= 5 and ("allele1" in low and "allele2" in low):
            return FormatType.ANCESTRYDNA, "tab-oddelena hlavicka s allele1/allele2 (AncestryDNA)"
        if is_header and "\t" in line and "rsid" in low and ("genotype" in low or "position" in low):
            return FormatType.TWENTYTHREEANDME, "tab-oddelena hlavicka rsid/chromosome/position/genotype (23andMe)"
        if not is_header and len(tab) == 5 and (tab[0].startswith("rs") or tab[0].startswith("i")):
            return FormatType.ANCESTRYDNA, "5 tab-oddelenych poli, rsID/iID v 1. stlpci (AncestryDNA-like)"
        if not is_header and len(tab) == 4 and (tab[0].startswith("rs") or tab[0].startswith("i")):
            return FormatType.TWENTYTHREEANDME, "4 tab-oddelene polia, rsID/iID v 1. stlpci (23andMe-like)"

    # --- 4. komentarove klucove slova (posledna moznost) ---
    if "ancestrydna" in lower or "ancestry.com" in lower:
        return FormatType.ANCESTRYDNA, "komentar obsahuje 'AncestryDNA'"
    if "familytreedna" in lower or "ftdna" in lower:
        return FormatType.FTDNA, "komentar obsahuje 'FamilyTreeDNA'"
    if "myheritage" in lower:
        return FormatType.MYHERITAGE, "komentar obsahuje 'MyHeritage'"
    if "23andme" in lower:
        return FormatType.TWENTYTHREEANDME, "komentar obsahuje '23andMe'"

    return FormatType.UNKNOWN, "nepodarilo sa rozpoznat ziaden znamy podpis formatu"


def detect_format(path: str) -> FormatType:
    return detect_format_verbose(path)[0]
