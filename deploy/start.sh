#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

export PUBLIC_MODE="${PUBLIC_MODE:-1}"
export REFERENCE_DIR="${REFERENCE_DIR:-/app/ref}"
export BEAGLE_JAR="${BEAGLE_JAR:-$REFERENCE_DIR/beagle.jar}"
export REF_DIR="${REF_DIR:-$REFERENCE_DIR/b37.bref3}"
export MAP_DIR="${MAP_DIR:-$REFERENCE_DIR/maps}"
export FASTA="${FASTA:-$REFERENCE_DIR/human_g1k_v37.fasta}"
export JOBS="${JOBS:-1}"
export XMX_JOB="${XMX_JOB:-3g}"

reference_ready() {
  local bref_count map_count
  bref_count="$(find "$REF_DIR" -maxdepth 1 -type f -name '*.bref3' 2>/dev/null | wc -l | tr -d '[:space:]')"
  map_count="$(find "$MAP_DIR" -maxdepth 1 -type f -name 'plink.chr*.GRCh37.map' 2>/dev/null | wc -l | tr -d '[:space:]')"
  [ -s "$BEAGLE_JAR" ] && [ -s "$FASTA" ] && [ -s "${FASTA}.fai" ] && \
    [ "${bref_count:-0}" -ge 23 ] && [ "${map_count:-0}" -ge 23 ]
}

if [ "${BOOTSTRAP_REFERENCE:-0}" = "1" ]; then
  if reference_ready; then
    echo "Railway: referenčný balík je pripravený v $REFERENCE_DIR"
  else
    echo "Railway: prvé spustenie — sťahujem a overujem referenčný balík"
    bash "$ROOT/setup_reference.sh"
  fi
fi

exec uvicorn app.main:app \
  --host 0.0.0.0 \
  --port "${PORT:-8000}" \
  --proxy-headers \
  --forwarded-allow-ips='*'
