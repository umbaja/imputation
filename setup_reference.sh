#!/usr/bin/env bash
# =============================================================================
# setup_reference.sh — stiahne VSETKO potrebne pre lokalnu Beagle imputaciu
# =============================================================================
# Stiahne: Beagle jar, geneticke mapy GRCh37, 1000G b37 bref3 panel (chr1-22),
# a referencnu FASTA GRCh37. Spusti RAZ na svojom stroji (ma pristup na internet).
#
# POUZITIE:  ./setup_reference.sh            # stiahne chr1-22
#            CHROMS="20 21 22" ./setup_reference.sh   # len vybrane chromozomy
# =============================================================================
set -euo pipefail
REF=ref
mkdir -p "$REF/b37.bref3" "$REF/maps"
CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"
BASE_BREF="https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.bref3"
BASE_MAP="https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps"

echo "[1/4] Beagle jar"
if [ ! -f "$REF/beagle.jar" ]; then
  # zisti aktualny nazov jar na stranke Beagle a stiahni; fallback na znamu verziu
  curl -fsSL -o "$REF/beagle.jar" \
    "https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar"
fi
echo "    -> $REF/beagle.jar"

echo "[2/4] Geneticke mapy GRCh37"
if [ ! -f "$REF/maps/plink.chr22.GRCh37.map" ]; then
  curl -fsSL -o "$REF/maps/plink.GRCh37.map.zip" "$BASE_MAP/plink.GRCh37.map.zip"
  unzip -o "$REF/maps/plink.GRCh37.map.zip" -d "$REF/maps" >/dev/null
fi
echo "    -> $REF/maps/plink.chrN.GRCh37.map"

echo "[3/4] 1000G b37 bref3 panel (chromozomy: $CHROMS)"
# Tip: chrX dotiahnes samostatne:  CHROMS="X" bash setup_reference.sh
dl_bref3() {                       # $1 = chromozom (1..22 alebo X)
  local c="$1"
  local f="$REF/b37.bref3/chr${c}.1kg.phase3.v5a.b37.bref3"
  [ -f "$f" ] && return 0
  echo "    sťahujem chr${c} ..."
  # X byva inak verzovany (v5a/v5b) — skus varianty, uloz vzdy pod v5a nazvom
  local name
  for name in "chr${c}.1kg.phase3.v5a.b37.bref3" \
              "chr${c}.1kg.phase3.v5b.b37.bref3"; do
    if curl -fsSL -o "$f" "$BASE_BREF/${name}"; then
      return 0
    fi
  done
  rm -f "$f"
  echo "    !! chr${c}: bref3 sa nepodarilo stiahnut. Skontroluj nazov v adresari:"
  echo "       $BASE_BREF/"
  echo "       a uloz ho ako $f"
  return 1
}
for c in $CHROMS; do
  dl_bref3 "$c" || echo "    (chr${c} preskocene)"
done
echo "    -> $REF/b37.bref3/"

echo "[4/4] Referencna FASTA GRCh37 (human_g1k_v37)"
if [ ! -f "$REF/human_g1k_v37.fasta" ]; then
  curl -fsSL -o "$REF/human_g1k_v37.fasta.gz" \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"
  gunzip -k "$REF/human_g1k_v37.fasta.gz"
  # index pre bcftools +fixref
  command -v samtools >/dev/null && samtools faidx "$REF/human_g1k_v37.fasta" || true
fi
echo "    -> $REF/human_g1k_v37.fasta"

echo "HOTOVO. Referencia je v ./$REF/"
echo "Nastav pred spustenim imputacie:"
echo "  export BEAGLE_JAR=\$PWD/$REF/beagle.jar"
echo "  export REF_DIR=\$PWD/$REF/b37.bref3"
echo "  export MAP_DIR=\$PWD/$REF/maps"
echo "  export FASTA=\$PWD/$REF/human_g1k_v37.fasta"
