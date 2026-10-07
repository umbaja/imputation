#!/usr/bin/env bash
# =============================================================================
# setup_wgs.sh — jednorazova priprava na spracovanie WGS/WES dat (GATK VCF + BAM)
# =============================================================================
# Doinstaluje pyliftover, stiahne UCSC chain subory (hg38<->hg19) a predpocita
# 23andMe sablonu v hg38 suradniciach (aby sa nemusela liftovat pri kazdom behu).
# Spustaj vo WSL z korena projektu:   bash setup_wgs.sh
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

VENV=.venv-wsl
CHAINS=ref/chains
IDREF=data/id_position_reference_v3v4v5.tsv.gz

echo ">> 1/4 virtualne prostredie"
if [ ! -f "$VENV/bin/activate" ]; then
  echo "   vytvaram $VENV ..."
  python3 -m venv "$VENV"
fi
# shellcheck disable=SC1090
source "$VENV/bin/activate"
pip install -q --upgrade pip
pip install -q pyliftover
python3 -c "import pyliftover; print('   pyliftover OK')"

echo ">> 2/4 chain subory z UCSC"
mkdir -p "$CHAINS"
base=https://hgdownload.soe.ucsc.edu/goldenPath
get() {  # url, cielovy subor
  if [ -s "$2" ]; then echo "   uz mam $(basename "$2")"; return; fi
  echo "   stahujem $(basename "$2") ..."
  curl -fSL --retry 3 -o "$2.part" "$1"
  mv "$2.part" "$2"
}
get "$base/hg38/liftOver/hg38ToHg19.over.chain.gz" "$CHAINS/hg38ToHg19.over.chain.gz"
get "$base/hg19/liftOver/hg19ToHg38.over.chain.gz" "$CHAINS/hg19ToHg38.over.chain.gz"
ls -lh "$CHAINS"

echo ">> 3/4 kontrola nastrojov"
if command -v samtools >/dev/null 2>&1; then
  echo "   samtools OK ($(samtools --version | head -1))"
else
  echo "   !! samtools CHYBA — bez neho sa neda citat pokrytie z BAM."
  echo "      Nainstaluj: sudo apt install -y samtools"
fi
if [ -f "${FASTA:-ref/human_g1k_v37.fasta}.fai" ]; then
  echo "   FASTA index OK"
else
  echo "   !! chyba ${FASTA:-ref/human_g1k_v37.fasta}.fai — sprav: samtools faidx ref/human_g1k_v37.fasta"
fi

echo ">> 4/4 predpocitanie 23andMe sablony v hg38 (raz, trva par minut)"
if [ ! -f "$IDREF" ]; then
  echo "   !! chyba $IDREF — sablonu neviem postavit."
  exit 1
fi
python3 - "$IDREF" <<'PY'
import sys, time
from dna_tools import wgs
t0 = time.time()
if wgs.TEMPLATE_CACHE.exists():
    print("   sablona uz existuje:", wgs.TEMPLATE_CACHE)
else:
    wgs.build_template_hg38(sys.argv[1], progress=lambda m: print("   ", m, flush=True))
    print("   hotovo za %.0f s" % (time.time() - t0))
print("   stav:", wgs.wgs_status())
PY

echo
echo ">> HOTOVO. Restartuj appku (pkill -9 -f uvicorn; bash start_app.sh)."
