"""
dna_tools
=========

Nástroje na detekciu formátu surových genetických dát, porovnanie zhodných
pozícií medzi rôznymi formátmi a prevod do univerzálneho 23andMe formátu.

Podporované vstupné formáty:
  - 23andMe            (raw data text export)
  - AncestryDNA        (raw data text export)
  - MyHeritage         (raw data CSV export)
  - FamilyTreeDNA      (raw data CSV export)
  - VCF                (Variant Call Format, v4.x)
  - PLINK              (.ped/.map textové súbory)
  - Illumina Final Report (GenomeStudio export)

Poznámka o vláknach (strand): 23andMe hlási genotypy na plus (forward) vlákne
referencie. Väčšina DTC formátov (23andMe, AncestryDNA, MyHeritage, FTDNA) už
je na plus vlákne. Illumina Final Report je typicky na TOP/BOT vlákne a
vyžaduje normalizáciu podľa manifestu čipu — konvertor na to upozorní.
"""

from .common import GenotypeData, NO_CALL
from .detect import detect_format, FormatType
from .parsers import parse_file
from .convert import to_23andme, write_23andme
from .compare import compare_files, ComparisonResult

__all__ = [
    "GenotypeData",
    "NO_CALL",
    "detect_format",
    "FormatType",
    "parse_file",
    "to_23andme",
    "write_23andme",
    "compare_files",
    "ComparisonResult",
]

__version__ = "0.1.0"
