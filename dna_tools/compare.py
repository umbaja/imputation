"""
compare.py
==========

Porovnanie zhodných pozícií a genotypov medzi dvoma (alebo viacerými) súbormi
v rôznych formátoch.

Čo počíta:
  - počet variantov v každom súbore,
  - prekryv pozícií podľa rsID a nezávisle podľa chr:pos,
  - konkordanciu genotypov na prekrývajúcich sa pozíciách:
        * exact      – identický genotyp (po zoradení alel, napr. AG == GA),
        * strand_flip– zhoda až po strand-flipe (napr. AG vs TC),
        * mismatch   – reálna nezhoda,
        * no_call    – aspoň jeden zo súborov má '--',
  - export nezhôd do CSV.

Použitie na "posúdenie zhodných pozícií v rôznych formátoch": porovná, do akej
miery ten istý vzorok exportovaný v dvoch formátoch pokrýva rovnaké SNP a či
sedia genotypy — kľúčové pri validácii konverzie a pred imputáciou.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import pandas as pd

from .common import (
    GenotypeData,
    NO_CALL,
    complement_genotype,
    sorted_genotype,
)
from .parsers import parse_file


@dataclass
class ComparisonResult:
    n_a: int
    n_b: int
    overlap_rsid: int
    overlap_pos: int
    exact: int
    strand_flip: int
    mismatch: int
    no_call: int
    merged: pd.DataFrame          # zlúčené prekrývajúce sa pozície + klasifikácia
    key: str                      # "rsid" alebo "pos"

    @property
    def concordance(self) -> float:
        """Podiel exact (+ strand_flip) zhôd z porovnateľných (ne-no-call) miest."""
        comparable = self.exact + self.strand_flip + self.mismatch
        if comparable == 0:
            return float("nan")
        return (self.exact + self.strand_flip) / comparable

    def summary(self) -> str:
        conc = self.concordance
        conc_s = "n/a" if conc != conc else f"{conc * 100:.2f}%"
        return (
            f"Súbor A: {self.n_a} variantov | Súbor B: {self.n_b} variantov\n"
            f"Prekryv podľa rsID: {self.overlap_rsid} | podľa chr:pos: {self.overlap_pos}\n"
            f"Porovnanie na kľúči '{self.key}':\n"
            f"  exact zhoda      : {self.exact}\n"
            f"  strand-flip zhoda: {self.strand_flip}\n"
            f"  nezhoda          : {self.mismatch}\n"
            f"  no-call (aspoň 1): {self.no_call}\n"
            f"  konkordancia     : {conc_s}"
        )


def _classify(gt_a: str, gt_b: str) -> str:
    if gt_a == NO_CALL or gt_b == NO_CALL:
        return "no_call"
    sa, sb = sorted_genotype(gt_a), sorted_genotype(gt_b)
    if sa == sb:
        return "exact"
    if sorted_genotype(complement_genotype(gt_a)) == sb:
        return "strand_flip"
    return "mismatch"


def compare(data_a: GenotypeData, data_b: GenotypeData, key: str = "rsid") -> ComparisonResult:
    """Porovná dve GenotypeData. key = 'rsid' alebo 'pos' (chr:pos)."""
    a = data_a.df.copy()
    b = data_b.df.copy()

    a["poskey"] = a["chromosome"].astype(str) + ":" + a["position"].astype("Int64").astype(str)
    b["poskey"] = b["chromosome"].astype(str) + ":" + b["position"].astype("Int64").astype(str)

    # prekryvy (na reporting)
    overlap_rsid = len(set(a["rsid"]) & set(b["rsid"]))
    overlap_pos = len(set(a["poskey"]) & set(b["poskey"]))

    join_col = "rsid" if key == "rsid" else "poskey"

    # deduplikuj podľa kľúča (ponechaj prvý výskyt)
    a_d = a.drop_duplicates(subset=join_col)
    b_d = b.drop_duplicates(subset=join_col)

    # stĺpce z A na ponechanie (bez duplikovania join_col)
    a_cols = [join_col] + [c for c in ["rsid", "chromosome", "position", "genotype"] if c != join_col]

    merged = pd.merge(
        a_d[a_cols],
        b_d[[join_col, "genotype"]],
        on=join_col, how="inner", suffixes=("_a", "_b"),
    )

    merged["class"] = [
        _classify(ga, gb) for ga, gb in zip(merged["genotype_a"], merged["genotype_b"])
    ]

    counts = merged["class"].value_counts().to_dict()

    return ComparisonResult(
        n_a=len(data_a.df),
        n_b=len(data_b.df),
        overlap_rsid=overlap_rsid,
        overlap_pos=overlap_pos,
        exact=counts.get("exact", 0),
        strand_flip=counts.get("strand_flip", 0),
        mismatch=counts.get("mismatch", 0),
        no_call=counts.get("no_call", 0),
        merged=merged,
        key=key,
    )


def compare_files(path_a: str, path_b: str, key: str = "rsid",
                  mismatch_csv: Optional[str] = None) -> ComparisonResult:
    """Načíta dva súbory (autodetekcia formátu) a porovná ich.

    mismatch_csv : ak je zadané, zapíše tam riadky klasifikované ako 'mismatch'.
    """
    data_a = parse_file(path_a)
    data_b = parse_file(path_b)
    result = compare(data_a, data_b, key=key)

    if mismatch_csv:
        mism = result.merged[result.merged["class"] == "mismatch"]
        mism.to_csv(mismatch_csv, index=False)

    return result
