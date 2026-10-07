#!/usr/bin/env bash
# start_app.sh — čistý (re)štart Genome Converter appky vo WSL.
# Zhodí starý uvicorn na porte 8000, doinštaluje závislosti a spustí nový server.
# Použitie (vo WSL, z koreňa projektu):   bash start_app.sh
set -e
cd "$(dirname "$0")"

PORT=8000

# 1) zhoď starý bežiaci server (aby sa zmeny v parsers.py/panel.py/main.py prejavili)
echo ">> Zhadzujem starý uvicorn na porte $PORT (ak beží)..."
pkill -f "uvicorn app.main:app" 2>/dev/null || true
if command -v fuser >/dev/null 2>&1; then
  fuser -k ${PORT}/tcp 2>/dev/null || true
fi
sleep 1

# 2) virtuálne prostredie + závislosti
# POZOR: .venv býva Windows venv (Scripts/, nie bin/) z run.bat/run.ps1.
# Pre WSL používame samostatný .venv-wsl, aby si prostredia nekolidovali.
VENV=.venv-wsl
if [ ! -f "$VENV/bin/activate" ]; then
  echo ">> Vytváram $VENV ..."
  rm -rf "$VENV"
  python3 -m venv "$VENV"
fi
source "$VENV/bin/activate"
pip install -q --upgrade pip
pip install -q -r requirements.txt

# 3) referencie pre imputáciu (krok 3 v UI) — bez nich ostáva tlačidlo šedé
export BEAGLE_JAR="$PWD/ref/beagle.jar"
export REF_DIR="$PWD/ref/b37.bref3"
export FASTA="$PWD/ref/human_g1k_v37.fasta"
# rýchla kontrola pripravenosti (informatívne)
[ -f "$BEAGLE_JAR" ] && echo ">> beagle.jar OK" || echo ">> CHÝBA $BEAGLE_JAR"
n_bref=$(ls "$REF_DIR"/*.bref3 2>/dev/null | wc -l); echo ">> bref3 súborov: $n_bref"
[ -f "$FASTA" ] && echo ">> fasta OK" || echo ">> CHÝBA $FASTA"
command -v java >/dev/null 2>&1 && echo ">> java OK ($(java -version 2>&1 | head -1))" \
  || echo ">> CHÝBA java — imputácia ostane vypnutá (nainštaluj: sudo apt install -y default-jre)"

# 4) spusti server
echo
echo ">> Appka beží na http://127.0.0.1:${PORT}   (Ctrl+C ukončí)"
echo ">> Illumina convert: zaškrtni 'online' aby sa 74 rs-only sond doplnilo cez Ensembl GRCh37."
exec python -m uvicorn app.main:app --host 127.0.0.1 --port ${PORT}
