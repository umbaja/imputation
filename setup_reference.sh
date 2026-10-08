#!/usr/bin/env bash
# Download and verify the external GRCh37/Beagle reference bundle.
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${REFERENCE_MANIFEST:-$ROOT/reference_manifest.tsv}"
REF="${REFERENCE_DIR:-$ROOT/ref}"
CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X}"
VERIFY_ONLY=0

usage() {
  cat <<'EOF'
Použitie:
  bash setup_reference.sh
  CHROMS="20 21 22" bash setup_reference.sh
  bash setup_reference.sh --verify-only

Premenné:
  CHROMS              chromozómy bref3, predvolene 1-22 a X
  REFERENCE_DIR       cieľový priečinok, predvolene ./ref
  REFERENCE_MANIFEST  iný manifest, predvolene ./reference_manifest.tsv
EOF
}

case "${1:-}" in
  "") ;;
  --verify-only) VERIFY_ONLY=1 ;;
  -h|--help) usage; exit 0 ;;
  *) usage >&2; exit 2 ;;
esac

for command_name in sha256sum wc tr awk; do
  command -v "$command_name" >/dev/null 2>&1 || {
    echo "CHYBA: chýba príkaz $command_name" >&2
    exit 1
  }
done
if [ "$VERIFY_ONLY" -eq 0 ]; then
  for command_name in curl unzip gzip; do
    command -v "$command_name" >/dev/null 2>&1 || {
      echo "CHYBA: chýba príkaz $command_name" >&2
      exit 1
    }
  done
fi
[ -f "$MANIFEST" ] || { echo "CHYBA: manifest neexistuje: $MANIFEST" >&2; exit 1; }
mkdir -p "$REF/b37.bref3" "$REF/maps"

verify_file() {
  local path="$1" expected_bytes="$2" expected_sha="$3"
  local actual_bytes actual_sha
  [ -f "$path" ] || return 1
  actual_bytes="$(wc -c < "$path" | tr -d '[:space:]')"
  [ "$actual_bytes" = "$expected_bytes" ] || return 1
  actual_sha="$(sha256sum "$path" | awk '{print $1}')"
  [ "$actual_sha" = "$expected_sha" ]
}

selected_chromosome() {
  local chromosome="$1" wanted
  for wanted in $CHROMS; do
    [ "$wanted" = "$chromosome" ] && return 0
  done
  return 1
}

download_verified() {
  local target="$1" expected_bytes="$2" expected_sha="$3" url="$4"
  local part backup
  if verify_file "$target" "$expected_bytes" "$expected_sha"; then
    echo "  OK (cache): ${target#"$ROOT/"}"
    return 0
  fi
  if [ "$VERIFY_ONLY" -eq 1 ]; then
    echo "  CHYBA: chýba alebo nesedí checksum: $target" >&2
    return 1
  fi

  mkdir -p "$(dirname -- "$target")"
  part="${target}.part"
  echo "  sťahujem: $url"
  curl --fail --location --retry 3 --retry-delay 2 --continue-at - \
    --output "$part" "$url"
  if ! verify_file "$part" "$expected_bytes" "$expected_sha"; then
    backup="${part}.invalid-$(date +%Y%m%d%H%M%S)-$$"
    mv -- "$part" "$backup"
    echo "  CHYBA: stiahnutý súbor nesedí s manifestom; zachovaný ako $backup" >&2
    return 1
  fi
  if [ -f "$target" ]; then
    backup="${target}.invalid-$(date +%Y%m%d%H%M%S)-$$"
    mv -- "$target" "$backup"
    echo "  pôvodný neplatný súbor: $backup"
  fi
  mv -- "$part" "$target"
  echo "  OK: ${target#"$ROOT/"}"
}

MAP_ARCHIVE=""
FASTA_ARCHIVE=""
FASTA_PATH=""
FASTA_BYTES=""
FASTA_SHA=""
BREF_COUNT=0

echo "Manifest: $MANIFEST"
echo "Referencia: $REF"
echo "Chromozómy: $CHROMS"
echo

while IFS=$'\t' read -r kind chromosome relative_path size_bytes sha256 url; do
  [ "$kind" = "kind" ] && continue
  [ -z "$kind" ] && continue
  case "$relative_path" in
    /*|../*|*/../*) echo "CHYBA: nebezpečná cesta v manifeste: $relative_path" >&2; exit 1 ;;
  esac
  target="$REF/$relative_path"
  case "$kind" in
    bref3)
      selected_chromosome "$chromosome" || continue
      download_verified "$target" "$size_bytes" "$sha256" "$url"
      BREF_COUNT=$((BREF_COUNT + 1))
      ;;
    beagle)
      download_verified "$target" "$size_bytes" "$sha256" "$url"
      ;;
    maps_archive)
      download_verified "$target" "$size_bytes" "$sha256" "$url"
      MAP_ARCHIVE="$target"
      ;;
    fasta_archive)
      download_verified "$target" "$size_bytes" "$sha256" "$url"
      FASTA_ARCHIVE="$target"
      ;;
    fasta)
      FASTA_PATH="$target"
      FASTA_BYTES="$size_bytes"
      FASTA_SHA="$sha256"
      ;;
    *) echo "CHYBA: neznámy typ v manifeste: $kind" >&2; exit 1 ;;
  esac
done < "$MANIFEST"

[ "$BREF_COUNT" -gt 0 ] || { echo "CHYBA: CHROMS nevybralo žiadny panel z manifestu" >&2; exit 1; }
[ -n "$MAP_ARCHIVE" ] || { echo "CHYBA: manifest nemá maps_archive" >&2; exit 1; }
[ -n "$FASTA_ARCHIVE" ] || { echo "CHYBA: manifest nemá fasta_archive" >&2; exit 1; }
[ -n "$FASTA_PATH" ] || { echo "CHYBA: manifest nemá rozbalenú fasta" >&2; exit 1; }

echo
echo "Genetické mapy"
if [ "$VERIFY_ONLY" -eq 0 ]; then
  unzip -oq "$MAP_ARCHIVE" -d "$REF/maps"
fi
for chromosome in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X; do
  [ -f "$REF/maps/plink.chr${chromosome}.GRCh37.map" ] || {
    echo "CHYBA: chýba rozbalená mapa chr${chromosome}" >&2
    exit 1
  }
done
echo "  OK: mapy GRCh37"

echo
echo "Referenčná FASTA"
if verify_file "$FASTA_PATH" "$FASTA_BYTES" "$FASTA_SHA"; then
  echo "  OK (cache): ${FASTA_PATH#"$ROOT/"}"
elif [ "$VERIFY_ONLY" -eq 1 ]; then
  echo "  CHYBA: chýba alebo nesedí checksum: $FASTA_PATH" >&2
  exit 1
else
  fasta_part="${FASTA_PATH}.part"
  if verify_file "$fasta_part" "$FASTA_BYTES" "$FASTA_SHA"; then
    echo "  OK (cache, rozbalené): ${fasta_part#"$ROOT/"}"
  else
    echo "  rozbaľujem: ${FASTA_ARCHIVE#"$ROOT/"}"
    gzip_status=0
    gzip -cd "$FASTA_ARCHIVE" > "$fasta_part" || gzip_status=$?
    if [ "$gzip_status" -ne 0 ]; then
      echo "  UPOZORNENIE: gzip skončil s kódom $gzip_status; výsledok musí prejsť kontrolou manifestu" >&2
    fi
  fi
  if ! verify_file "$fasta_part" "$FASTA_BYTES" "$FASTA_SHA"; then
    echo "CHYBA: rozbalená FASTA nesedí s manifestom: $fasta_part" >&2
    exit 1
  fi
  if [ -f "$FASTA_PATH" ]; then
    fasta_backup="${FASTA_PATH}.invalid-$(date +%Y%m%d%H%M%S)-$$"
    mv -- "$FASTA_PATH" "$fasta_backup"
    echo "  pôvodná neplatná FASTA: $fasta_backup"
  fi
  mv -- "$fasta_part" "$FASTA_PATH"
  echo "  OK: ${FASTA_PATH#"$ROOT/"}"
fi

if [ "$VERIFY_ONLY" -eq 0 ] && command -v samtools >/dev/null 2>&1; then
  if [ ! -f "${FASTA_PATH}.fai" ] || [ "$FASTA_PATH" -nt "${FASTA_PATH}.fai" ]; then
    samtools faidx "$FASTA_PATH"
  fi
fi

echo
if [ "$VERIFY_ONLY" -eq 1 ]; then
  echo "OVERENÉ. Vybraná referencia zodpovedá manifestu."
else
  echo "HOTOVO. Referencia bola stiahnutá a overená podľa manifestu."
fi
echo "Nastav pred spustením imputácie:"
echo "  export BEAGLE_JAR=\"$REF/beagle.jar\""
echo "  export REF_DIR=\"$REF/b37.bref3\""
echo "  export MAP_DIR=\"$REF/maps\""
echo "  export FASTA=\"$FASTA_PATH\""
