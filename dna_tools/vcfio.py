"""
vcfio.py — 23andMe -> VCF. Dva rezimy:
  A) --fasta: REF z FASTA, emituje len informativne (het/hom-alt) SNP.
  B) --sites: REF/ALT z panelovych "sites" (chrom,pos,ref,alt) -> emituje AJ
     homozygotne-ref markery (0/0), cim sa vyrazne zvysi pocet kotiev pre Beagle.
Plus extract_sites_vcf(): vytiahne (chrom,pos,ref,alt) z VCF (napr. z imputovaneho
vystupu 1. behu Beagle) -> panelove sites bez extra stahovania.
"""

from __future__ import annotations

import gzip
from typing import Dict, List, Optional, Tuple

from .parsers import parse_file

_ACGT = set("ACGT")


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_fasta_contig(fasta_path: str, chrom: str) -> str:
    want = {str(chrom), "chr" + str(chrom)}
    parts: List[str] = []
    capturing = False
    with open(fasta_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if capturing:
                    break
                capturing = line[1:].split()[0] in want
                continue
            if capturing:
                parts.append(line.strip())
    return "".join(parts).upper()


def extract_sites_vcf(vcf_path: str) -> Dict[Tuple[str, int], Tuple[str, str]]:
    """Z VCF vytiahne bialelicke SNP sites: {(chrom,pos): (ref,alt)}."""
    sites: Dict[Tuple[str, int], Tuple[str, str]] = {}
    with _open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, pos, ref, alt = f[0], f[1], f[3].upper(), f[4].upper()
            if len(ref) != 1 or ref not in _ACGT:
                continue
            if len(alt) != 1 or alt not in _ACGT:
                continue  # len bialelicke SNP
            try:
                sites[(str(chrom), int(pos))] = (ref, alt)
            except ValueError:
                pass
    return sites


def write_vcf_from_sample(sample_path, out_vcf, fasta_path=None, sites=None,
                          chroms=None, sample_name="SAMPLE"):
    """Zapise VCF z 23andMe vzorky.

    Rezim A (fasta_path): REF z FASTA, len informativne SNP.
    Rezim B (sites): dict {(chrom,pos):(ref,alt)} -> emituje aj hom-ref (viac kotiev).
    chroms obmedzi na dane chromozomy (list/str); ak None, vsetky.
    """
    data = parse_file(sample_path)
    df = data.df
    if chroms is not None:
        if isinstance(chroms, (str, int)):
            chroms = [str(chroms)]
        chroms = set(str(c) for c in chroms)
    stats = {"written": 0, "skipped_homref": 0, "skipped_nonsnp": 0,
             "skipped_mismatch": 0, "skipped_nosite": 0}
    lines = []
    contigs = set()

    fasta_cache: Dict[str, str] = {}

    for row in df.itertuples(index=False):
        chrom = str(row.chromosome)
        if chroms is not None and chrom not in chroms:
            continue
        gt = str(row.genotype)
        try:
            pos = int(row.position)
        except (TypeError, ValueError):
            continue
        if len(gt) != 2 or gt[0] not in _ACGT or gt[1] not in _ACGT:
            stats["skipped_nonsnp"] += 1
            continue
        a1, a2 = gt[0], gt[1]

        if sites is not None:
            key = (chrom, pos)
            site = sites.get(key)
            if site is None:
                stats["skipped_nosite"] += 1
                continue
            ref, alt = site
            valid = {ref, alt}
            if a1 not in valid or a2 not in valid:
                stats["skipped_mismatch"] += 1   # ina alela / strand
                continue
        else:
            seq = fasta_cache.get(chrom)
            if seq is None:
                seq = read_fasta_contig(fasta_path, chrom)
                fasta_cache[chrom] = seq
            if not seq or pos < 1 or pos > len(seq):
                continue
            ref = seq[pos - 1]
            if ref not in _ACGT:
                stats["skipped_mismatch"] += 1
                continue
            non_ref = sorted(a for a in {a1, a2} if a != ref)
            if len(non_ref) == 0:
                stats["skipped_homref"] += 1
                continue
            if len(non_ref) == 2:
                stats["skipped_mismatch"] += 1
                continue
            alt = non_ref[0]

        idx = {ref: 0, alt: 1}
        g = "%d/%d" % (idx[a1], idx[a2])
        rsid = str(row.rsid) if str(row.rsid).startswith("rs") else "."
        lines.append((chrom, pos, "%s\t%d\t%s\t%s\t%s\t.\t.\t.\tGT\t%s"
                      % (chrom, pos, rsid, ref, alt, g)))
        contigs.add(chrom)
        stats["written"] += 1

    def _ck(c):
        try:
            return (0, int(c))
        except ValueError:
            return (1, c)
    lines.sort(key=lambda t: (_ck(t[0]), t[1]))

    with open(out_vcf, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c in sorted(contigs, key=_ck):
            fh.write("##contig=<ID=%s>\n" % c)
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % sample_name)
        for _, _, ln in lines:
            fh.write(ln + "\n")
    return stats
