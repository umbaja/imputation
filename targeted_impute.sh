#!/usr/bin/env bash
# =============================================================================
# targeted_impute.sh — CIELENA imputacia len panelovych lokusov (nie cely cip)
# =============================================================================
# Okolo cielov (impute_targets.csv) vyrezze uzke regiony a imputuje LEN tie
# pomocou Beagle 5 + 1000G. Bezi LOKALNE. Paralelne po regionoch + durable cache.
# POUZITIE:
#   ./targeted_impute.sh SAMPLE_23ANDME.txt TARGETS.csv PANEL_STATUS.csv OUT_PREFIX
# =============================================================================
set -euo pipefail

export LC_ALL=C
export LANG=C

SAMPLE="${1:?vstupny normalizovany genotypovy subor}"
TARGETS="${2:?CSV s cielmi (rsid,chromosome,position)}"
STATUS="${3:?panel-check --out-all CSV}"
OUT="${4:-targeted}"

BEAGLE_JAR="${BEAGLE_JAR:-beagle.jar}"
REF_DIR="${REF_DIR:-ref/b37.bref3}"
MAP_DIR="${MAP_DIR:-ref/maps}"
FASTA="${FASTA:-ref/human_g1k_v37.fasta}"
FLANK="${FLANK:-250000}"
DR2MIN="${DR2MIN:-0.3}"
XMX="${XMX:-8g}"
# Timeout na JEDEN region (s). 0 = vypnute. Zaseknuty region sa zabije a preskoci.
REGION_TIMEOUT="${REGION_TIMEOUT:-0}"
# Paralelizmus: kolko Beagle behov naraz. Default min(nproc-2, 6).
_NPROC=$(nproc 2>/dev/null || echo 4)
if [ "$_NPROC" -gt 8 ]; then
  _JOBS_DEF=6
elif [ "$_NPROC" -gt 3 ]; then
  _JOBS_DEF=$((_NPROC - 2))
else
  _JOBS_DEF=1
fi
JOBS="${JOBS:-$_JOBS_DEF}"
XMX_JOB="${XMX_JOB:-3g}"
WORK_DIR="${WORK_DIR:-work_ti}"

mkdir -p "$WORK_DIR" "${OUT}_regions"

echo "[1/6] Cielova vzorka -> VCF (bez plink, cez dna_tools/vcfio; REF z FASTA)"
python3 - "$SAMPLE" "$FASTA" "$TARGETS" "$WORK_DIR" <<'PY'
import sys, csv
from pathlib import Path
from dna_tools.vcfio import write_vcf_from_sample
sample, fasta, targets, work_dir = sys.argv[1:5]
chroms = sorted({r["chromosome"] for r in csv.DictReader(open(targets)) if r.get("chromosome")})
stats = write_vcf_from_sample(sample, str(Path(work_dir) / "raw.vcf"), fasta_path=fasta,
                              chroms=chroms, sample_name="sample")
print("    chromozomy cielov:", ",".join(chroms))
print("    VCF:", stats)
PY
bgzip -f "$WORK_DIR/raw.vcf"

echo "[2/6] Sort + index (REF uz zarovnany na FASTA)"
bcftools sort "$WORK_DIR/raw.vcf.gz" -Oz -o "$WORK_DIR/sorted.vcf.gz"
bcftools index -t "$WORK_DIR/sorted.vcf.gz"

echo "[kontrola] duplicity cielovych lokusov (dedup podla: ${DEDUP_BY:-position})"
python3 - "$TARGETS" "${DEDUP_BY:-position}" "$WORK_DIR/targets_clean.csv" <<'PY'
import sys, csv
from collections import Counter, defaultdict
inp, mode, outp = sys.argv[1], sys.argv[2], sys.argv[3]
rows = [r for r in csv.DictReader(open(inp)) if r.get("position") not in (None, "")]
def poskey(r): return (str(r.get("chromosome", "")).strip(), str(int(float(r["position"]))))
def rskey(r):  return (r.get("rsid") or r.get("id") or "").strip()
n = len(rows)
pos_c = Counter(poskey(r) for r in rows)
rs_c  = Counter(rskey(r) for r in rows if rskey(r))
dup_pos = sum(v - 1 for v in pos_c.values() if v > 1)
dup_rs  = sum(v - 1 for v in rs_c.values() if v > 1)
key = poskey if mode == "position" else rskey
seen, out = set(), []
for r in rows:
    k = key(r)
    if k in seen: continue
    seen.add(k); out.append(r)
fields = list(rows[0].keys()) if rows else ["rsid", "chromosome", "position"]
with open(outp, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields); w.writeheader()
    for r in out: w.writerow(r)
print(f"    vstupnych cielov: {n}")
print(f"    duplicitne pozicie: {dup_pos} | duplicitne rsID: {dup_rs}")
print(f"    po deduplikacii ({mode}) zostava: {len(out)} cielov")
PY
TARGETS="$WORK_DIR/targets_clean.csv"

echo "[3/6] Vypocet regionov (+-${FLANK} bp) okolo cielov"
python3 - "$TARGETS" "$FLANK" > "$WORK_DIR/regions.tsv" <<'PY'
import sys, csv
from dna_tools.targeted_impute import build_windows
targets=[(r["chromosome"], int(float(r["position"])))
         for r in csv.DictReader(open(sys.argv[1]))
         if r.get("position") not in (None, "")]
for reg in build_windows(targets, flank=int(sys.argv[2])):
    print(f"{reg.chrom}\t{reg.start}\t{reg.end}")
PY
echo "    regionov: $(wc -l < "$WORK_DIR/regions.tsv")"

if [ ! -s "$WORK_DIR/regions.tsv" ]; then
  echo "CHYBA: 0 cielovych regionov — ziaden ciel nema pouzitelnu poziciu." >&2
  exit 3
fi

# --- CACHE ---
CACHE_DIR="${CACHE_DIR:-work/impute_cache}"
mkdir -p "$CACHE_DIR"
CACHE_KEY=$( { bcftools view -H "$WORK_DIR/sorted.vcf.gz"; cat "$WORK_DIR/regions.tsv"; } \
             | sha256sum | cut -c1-40 )
CACHED="$CACHE_DIR/${CACHE_KEY}_imputed.vcf.gz"
echo "    cache kluc: $CACHE_KEY"

# Durable per-region cache: kluc = suradnice regionu + GENOTYPY V TOM REGIONE
# (presne to, co Beagle dostane na vstup). Odolne voci kozmetickym/nedeterministickym
# rozdielom inde vo vzorke (napr. online Ensembl anotacia pri konverzii) — hitne
# vzdy, ak je vstup regionu rovnaky, aj ked sa zvysok suboru lisi.
RCACHE="${CACHE_DIR}/regions"
mkdir -p "$RCACHE"

if [ -f "$CACHED" ]; then
  echo "[4/6] CACHE HIT — tato vzorka+ciele uz boli imputovane, preskakujem Beagle"
  cp "$CACHED" "${OUT}_imputed.vcf.gz"
  bcftools index -f -t "${OUT}_imputed.vcf.gz"
  echo "[5/6] (preskocene — pouzity ulozeny vysledok)"
else
  echo "[4/6] CACHE MISS — spustam Beagle imputaciu po regionoch (paralelne: JOBS=${JOBS}, Xmx=${XMX_JOB}/job)"
  : > "$WORK_DIR/imputed_list.txt"

  impute_one() {
    local chr="$1" start="$2" end="$3"
    local refp="${REF_DIR}/chr${chr}.1kg.phase3.v5a.b37.bref3"
    local mapp="${MAP_DIR}/plink.chr${chr}.GRCh37.map"
    local out="${OUT}_regions/reg_${chr}_${start}_${end}"
    if [ ! -f "$refp" ] || [ ! -f "$mapp" ]; then
      echo "    preskakujem chr${chr}:${start}-${end} — chyba 1000G referencia/mapa (mimo chr1-22?)"
      echo ">> REGION_DONE chr${chr}:${start}-${end} skip"
      return 0
    fi
    if [ -s "${out}.vcf.gz" ] && bcftools index -f -t "${out}.vcf.gz" 2>/dev/null; then
      echo "    región chr${chr}:${start}-${end} uz hotovy — pouzijem (resume)"
      echo ">> REGION_DONE chr${chr}:${start}-${end} resume"
      return 0
    fi
    local rkey rcf
    rkey=$( { printf '%s:%s:%s\n' "$chr" "$start" "$end"; \
              bcftools view -H -r "${chr}:${start}-${end}" "$WORK_DIR/sorted.vcf.gz" 2>/dev/null; } \
            | sha256sum | cut -c1-40)
    rcf="${RCACHE}/${rkey}.vcf.gz"
    if [ -s "$rcf" ]; then
      cp "$rcf" "${out}.vcf.gz"
      bcftools index -f -t "${out}.vcf.gz" 2>/dev/null || true
      echo "    región chr${chr}:${start}-${end} z cache (uz raz zbehol)"
      echo ">> REGION_DONE chr${chr}:${start}-${end} cache"
      return 0
    fi
    local rc_run=0
    if [ "${REGION_TIMEOUT}" != "0" ] && command -v timeout >/dev/null 2>&1; then
      timeout "${REGION_TIMEOUT}" java -Xmx"${XMX_JOB}" -jar "$BEAGLE_JAR" \
        gt="$WORK_DIR/sorted.vcf.gz" ref="$refp" map="$mapp" \
        chrom="${chr}:${start}-${end}" impute=true gp=true \
        out="${out}" >/dev/null 2>&1 || rc_run=$?
    else
      java -Xmx"${XMX_JOB}" -jar "$BEAGLE_JAR" \
        gt="$WORK_DIR/sorted.vcf.gz" ref="$refp" map="$mapp" \
        chrom="${chr}:${start}-${end}" impute=true gp=true \
        out="${out}" >/dev/null 2>&1 || rc_run=$?
    fi
    if [ "$rc_run" -ne 0 ]; then
      echo "    Beagle zlyhal/timeout (rc=$rc_run) na chr${chr}:${start}-${end}, preskakujem"
      rm -f "${out}.vcf.gz" "${out}.vcf.gz.tbi"
      echo ">> REGION_DONE chr${chr}:${start}-${end} fail"
      return 0
    fi
    bcftools index -f -t "${out}.vcf.gz" 2>/dev/null || true
    cp "${out}.vcf.gz" "$rcf" 2>/dev/null || true
    echo "    región chr${chr}:${start}-${end} OK"
    echo ">> REGION_DONE chr${chr}:${start}-${end} ok"
    return 0
  }

  while IFS=$'\t' read -r chr start end; do
    [ -z "${chr:-}" ] && continue
    impute_one "$chr" "$start" "$end" &
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do
      wait -n 2>/dev/null || sleep 0.2
    done
  done < "$WORK_DIR/regions.tsv"
  wait

  total=$(wc -l < "$WORK_DIR/regions.tsv")
  while IFS=$'\t' read -r chr start end; do
    [ -z "${chr:-}" ] && continue
    out="${OUT}_regions/reg_${chr}_${start}_${end}"
    if [ -s "${out}.vcf.gz" ] && [ -f "${out}.vcf.gz.tbi" ]; then
      echo "${out}.vcf.gz" >> "$WORK_DIR/imputed_list.txt"
    fi
  done < "$WORK_DIR/regions.tsv"
  done_n=$(wc -l < "$WORK_DIR/imputed_list.txt")
  echo "    preskocenych regionov: $(( total - done_n ))"

  if [ ! -s "$WORK_DIR/imputed_list.txt" ]; then
    echo "CHYBA: ani jeden region sa nenaimputoval (vsetky preskocene?)." >&2
    exit 4
  fi

  echo "[5/6] Spojenie regionov"
  bcftools concat -f "$WORK_DIR/imputed_list.txt" -a -Oz -o "${OUT}_imputed.vcf.gz"
  bcftools index -t "${OUT}_imputed.vcf.gz"
  cp "${OUT}_imputed.vcf.gz" "$CACHED"
  echo "    ulozene do cache: $CACHED"
fi

echo "[6/6] Extrakcia cielov (GT+DR2) a zlucenie do finalneho panela"
python3 - "$TARGETS" "$STATUS" "${OUT}_imputed.vcf.gz" "$DR2MIN" "${OUT}_panel_final.csv" <<'PY'
import sys, csv
from dna_tools.targeted_impute import extract_targets_from_vcf, merge_panel
targets_csv, status_csv, vcf, dr2min, outcsv = sys.argv[1:6]
targets={(r["chromosome"], int(float(r["position"]))): r["rsid"]
         for r in csv.DictReader(open(targets_csv))
         if r.get("position") not in (None, "")}
imp=extract_targets_from_vcf(vcf, targets, min_gp=float(dr2min))
final=merge_panel(status_csv, imputed=imp)
final["position"]=final["position"].astype(float).astype("Int64")
final.to_csv(outcsv, index=False)
vc=final["source"].value_counts().to_dict()
print("  zdroje:", vc)
print("  nizka spolahlivost (max_gp<%s):"%dr2min, int(imp["low_conf"].sum()))
print("  -> "+outcsv)

import pandas as _pd
miss=final[final["source"]=="missing"].copy()
def _cat(rs, ch, pos):
    rs=str(rs); ch=str(ch)
    if pos is None or (isinstance(pos,float) and _pd.isna(pos)) or str(pos) in ("","nan","<NA>"):
        return ("bez_pozicie", "nepodarilo sa urcit poziciu")
    if ch in ("X","Y","MT"):
        return ("chr_bez_referencie", "dotiahni 1000G referenciu pre chr%s"%ch)
    if rs.lower().startswith("i"):
        return ("23andMe_custom_proba", "firemna proba mimo 1000G")
    return ("mimo_1000G", "variant nie je v 1000G paneli")
unimp=outcsv.replace("_panel_final.csv","_unimputed.csv")
cols=["rsid","chromosome","position"]
if len(miss):
    cats=[_cat(r.rsid, r.chromosome, r.position) for r in miss.itertuples(index=False)]
    miss["kategoria"]=[c[0] for c in cats]
    miss["poznamka"]=[c[1] for c in cats]
    miss[cols+["kategoria","poznamka"]].to_csv(unimp, index=False)
    from collections import Counter
    print("  neimputovanych:", len(miss), dict(Counter(miss["kategoria"])), "->", unimp)
else:
    _pd.DataFrame(columns=cols+["kategoria","poznamka"]).to_csv(unimp, index=False)
    print("  neimputovanych: 0 ->", unimp)
PY

echo "HOTOVO. Finalny panel: ${OUT}_panel_final.csv"
