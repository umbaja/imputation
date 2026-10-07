"""
common.py
=========

Spolocna interna reprezentacia genotypovych dat a pomocne funkcie na
normalizaciu chromozomov, genotypov a alel.

Interna reprezentacia je pandas DataFrame so standardizovanymi stlpcami:
    rsid, chromosome, position, genotype
genotype: 1-2 znaky; diploid "AG", haploid "A", indel "II"/"D"; no-call = "--".
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import pandas as pd

NO_CALL = "--"
CANONICAL_CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "-": "-", "0": "0", "I": "I", "D": "D"}
COLUMNS = ["rsid", "chromosome", "position", "genotype"]


def normalize_chromosome(value) -> str:
    """Zjednoti zapis chromozomu na kanonicku formu (chr1->1, 23->X, 26->MT...)."""
    if value is None:
        return "0"
    s = str(value).strip()
    if s == "":
        return "0"
    if s.lower().startswith("chr"):
        s = s[3:]
    s = s.upper()
    mapping = {"23": "X", "24": "Y", "25": "X", "26": "MT", "M": "MT", "MITO": "MT", "XY": "X"}
    return mapping.get(s, s)


def _clean_allele(a) -> str:
    """Vycisti jednu alelu; neplatne/chybajuce -> prazdny retazec. Zachova A/C/G/T aj I/D."""
    if a is None:
        return ""
    a = str(a).strip().upper()
    if a in ("", "0", "-", ".", "N"):
        return ""
    return a[0]


def normalize_genotype(a1, a2=None) -> str:
    """Znormalizuje genotyp: no-call -> '--', diploid -> 2 znaky, haploid -> 1 znak."""
    if a2 is None:
        s = "" if a1 is None else str(a1).strip().upper()
        if s in ("", "--", "..", "00", "NN", "N", "-", "."):
            return NO_CALL
        if len(s) == 1:
            c = _clean_allele(s)
            return c if c else NO_CALL
        if len(s) == 2:
            b1, b2 = _clean_allele(s[0]), _clean_allele(s[1])
            if not b1 and not b2:
                return NO_CALL
            if b1 and b2:
                return b1 + b2
            return b1 or b2
        cleaned = "".join(_clean_allele(ch) for ch in s)
        return cleaned[:2] if cleaned else NO_CALL
    b1, b2 = _clean_allele(a1), _clean_allele(a2)
    if not b1 and not b2:
        return NO_CALL
    if b1 and b2:
        return b1 + b2
    return b1 or b2


def complement_genotype(genotype: str) -> str:
    """Komplement genotypu (strand-flip). 'AG' -> 'TC'."""
    return "".join(COMPLEMENT.get(b, b) for b in genotype)


def sorted_genotype(genotype: str) -> str:
    """Alely v abecednom poradi ('AG'=='GA'). No-call a haploid ostavaju nezmenene."""
    if genotype == NO_CALL or len(genotype) != 2:
        return genotype
    return "".join(sorted(genotype))


def chrom_sort_key(chrom: str):
    """Kluc na zoradenie chromozomov 1..22, X, Y, MT."""
    try:
        return (0, int(chrom))
    except (ValueError, TypeError):
        order = {"X": 23, "Y": 24, "MT": 25}
        return (0, order.get(chrom, 99))


@dataclass
class GenotypeData:
    """Kontajner pre parsovane genotypove data plus metadata."""

    df: pd.DataFrame
    source_format: str = "unknown"
    build: str = "GRCh37"
    sample_id: Optional[str] = None
    warnings: List[str] = field(default_factory=list)

    def __post_init__(self):
        missing = [c for c in COLUMNS if c not in self.df.columns]
        if missing:
            raise ValueError("Chybajuce stlpce v internej reprezentacii: %s" % missing)
        self.df = self.df[COLUMNS].copy()
        self.df["position"] = pd.to_numeric(self.df["position"], errors="coerce").astype("Int64")

    def sort(self) -> "GenotypeData":
        d = self.df.copy()
        d["_ck"] = d["chromosome"].map(chrom_sort_key)
        d = d.sort_values(["_ck", "position"]).drop(columns="_ck").reset_index(drop=True)
        self.df = d
        return self

    @property
    def n_variants(self) -> int:
        return len(self.df)

    @property
    def n_called(self) -> int:
        return int((self.df["genotype"] != NO_CALL).sum())

    def summary(self) -> str:
        return (
            "Format: %s | build: %s | vzorka: %s | variantov: %d | zavolanych: %d | no-call: %d"
            % (self.source_format, self.build, self.sample_id or "n/a",
               self.n_variants, self.n_called, self.n_variants - self.n_called)
        )
