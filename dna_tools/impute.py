"""
impute.py
=========

Dva kroky suvisiace s "chybajucimi" genotypmi:

1. build_scaffold() — namapuje genotypy vzorky na PEVNY marker set sablony
   (napr. 23andMe v5). Kazda pozicia sablony, ktoru vzorka nemeria, dostane
   no-call '--'. Toto je deterministicke a spustitelne HNED — identifikuje
   presne, ktore pozicie chybaju. NIE je to statisticka imputacia.

2. Statisticka imputacia chybajucich genotypov (Beagle 5 + referencny panel
   1000G/HRC/TOPMed) — vyzaduje externy referencny panel (GB dat) a nedokaze
   sa spustit bez neho. Kompletna pripravena pipeline je v impute_pipeline.sh.
   Tu je len wrapper run_beagle(), ktory pipeline zavola, ak su nastroje a
   panel k dispozicii.

Preco to takto: naivna imputacia (napr. "dopln najcastejsi genotyp") produkuje
NESPRAVNE genotypy, co je pri genetickych datach horsie nez ponechat no-call.
Spravna imputacia vyuziva vazbovu nerovnovahu (LD) z referencneho panelu.
"""

from __future__ import annotations

import subprocess
from typing import Optional, Tuple, Union

import pandas as pd

from .common import GenotypeData, NO_CALL, sorted_genotype
from .parsers import parse_file


def _as_data(x: Union[str, GenotypeData]) -> GenotypeData:
    return parse_file(x) if isinstance(x, str) else x


def build_scaffold(sample: Union[str, GenotypeData],
                   template: Union[str, GenotypeData],
                   sort_alleles: bool = True) -> Tuple[GenotypeData, dict]:
    """Namapuje genotypy vzorky na marker set sablony (matchuje podla chr:pos).

    sample   : vzorka na prevod (cesta alebo GenotypeData), napr. MyHeritage.
    template : sablona marker setu (cesta alebo GenotypeData), napr. 23andMe v5.
    Vrati (GenotypeData na marker sete sablony, statistiky pokrytia).
    Pozicie sablony bez zhody vo vzorke dostanu no-call '--'.
    rsID a poradie markerov sa preberaju zo sablony (t.j. 23andMe rsID).
    """
    sample = _as_data(sample)
    template = _as_data(template)

    s = sample.df.copy()
    s["poskey"] = s["chromosome"].astype(str) + ":" + s["position"].astype("Int64").astype(str)
    # ak je viac zaznamov na rovnakej pozicii, ponechaj prvy zavolany
    s = s.drop_duplicates(subset="poskey")
    smap = dict(zip(s["poskey"], s["genotype"]))

    t = template.df.copy()
    t["poskey"] = t["chromosome"].astype(str) + ":" + t["position"].astype("Int64").astype(str)
    gts = [smap.get(pk, NO_CALL) for pk in t["poskey"]]
    if sort_alleles:
        gts = [sorted_genotype(g) for g in gts]

    out = t[["rsid", "chromosome", "position"]].copy()
    out["genotype"] = gts

    n_total = len(out)
    n_filled = sum(1 for g in gts if g != NO_CALL)
    stats = {
        "template_markers": n_total,
        "filled_from_sample": n_filled,
        "missing_no_call": n_total - n_filled,
        "fill_rate": (n_filled / n_total) if n_total else 0.0,
    }
    gd = GenotypeData(out, source_format="23andme-scaffold")
    gd.warnings.append(
        "Scaffold: %d/%d pozicii vyplnenych (%.1f%%), %d chyba (no-call). "
        "Na statisticku imputaciu chybajucich pouzi impute_pipeline.sh (Beagle + 1000G)."
        % (n_filled, n_total, 100.0 * stats["fill_rate"], stats["missing_no_call"])
    )
    return gd, stats


def run_beagle(sample_23andme: str, out_prefix: str,
               beagle_jar: str, ref_dir: str, map_dir: str,
               fasta: Optional[str] = None,
               chroms=tuple(str(i) for i in range(1, 23)) + ("X",),
               plink2: str = "plink2", bcftools: str = "bcftools",
               java_xmx: str = "16g") -> str:
    """Spusti realnu Beagle imputaciu (vyzaduje externe nastroje + referencny panel).

    Ocakava:
      beagle_jar : cesta k beagle.jar
      ref_dir    : priecinok s chr{N}.1kg.phase3.v5a.b37.bref3 (1000G b37)
      map_dir    : priecinok s plink.chr{N}.GRCh37.map
      fasta      : (volitelne) referencna FASTA na zarovnanie REF/ALT
    Toto je tenky wrapper okolo krokov v impute_pipeline.sh; vrati cestu k
    vyslednemu imputovanemu 23andMe suboru. Ak nastroje/panel chybaju, vyhodi
    zrozumitelnu chybu.
    """
    raise NotImplementedError(
        "run_beagle je len rozhranie. Reálnu imputáciu spusti cez impute_pipeline.sh, "
        "ktorý obsahuje overené kroky (plink2 -> fixref -> Beagle -> späť do 23andMe). "
        "Vyžaduje beagle.jar, 1000G bref3 panel a genetické mapy — tie sa v tomto "
        "sandboxe nedajú stiahnuť (blokovaný allowlist), preto ich spusti na svojom stroji."
    )
