#!/usr/bin/env bash
# Stiahne VSETKY referencie pre plnu imputaciu: 1000G bref3 chr1-22 + FASTA chr1-22.
# Robustne: -c (resume), timeout, vela pokusov, preskoci hotove (.done marker).
# Spusti a nechaj bezat (aj hodinu). Progres sa loguje.
set -u
cd "$(dirname "$0")"
mkdir -p ref/b37.bref3 logs
LOG=logs/download_all.log
exec > >(tee "$LOG") 2>&1
CHROMS="${CHROMS:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"
W="wget -q -c --timeout=60 --tries=30 --retry-connrefused --waitretry=8"

echo "=== $(date '+%H:%M:%S') START (chromozomy: $CHROMS) ==="
[ -f ref/beagle.jar.done ] || { $W https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar -O ref/beagle.jar && touch ref/beagle.jar.done; }

for c in $CHROMS; do
  BREF=ref/b37.bref3/chr${c}.1kg.phase3.v5a.b37.bref3
  if [ -f "${BREF}.done" ]; then echo "chr$c bref3: HOTOVE (skip)"; else
    echo "chr$c bref3: stahujem ($(date '+%H:%M:%S')) ..."
    $W "http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.bref3/chr${c}.1kg.phase3.v5a.b37.bref3" -O "$BREF" \
      && touch "${BREF}.done" && echo "  chr$c bref3 OK: $(du -h $BREF|cut -f1)" || echo "  chr$c bref3 ZLYHAL (skus znova)"
  fi
  FA=ref/chr${c}.fa
  if [ -f "${FA}.done" ]; then echo "chr$c fasta: HOTOVE (skip)"; else
    echo "chr$c fasta: stahujem ..."
    if $W "https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.${c}.fa.gz" -O ${FA}.gz; then
      gunzip -f ${FA}.gz && touch "${FA}.done" && echo "  chr$c fasta OK: $(du -h $FA|cut -f1)"
    else echo "  chr$c fasta ZLYHAL"; fi
  fi
done
echo "=== $(date '+%H:%M:%S') HOTOVO. Skontroluj, ci su vsetky .done: ==="
ls ref/b37.bref3/*.done 2>/dev/null | wc -l | xargs echo "bref3 hotovych:"
ls ref/*.fa.done 2>/dev/null | wc -l | xargs echo "fasta hotovych:"
