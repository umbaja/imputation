"""
batch.py — DAVKOVE spracovanie: jedno tlacidlo spusti celu pipeline pre N vzoriek.

Pre kazdu vstupnu vzorku bezi:
    1) prevod do normalizovaneho genotypu (vratane Illumina relabel + rs-only anotacie)
    2) kontrola pokrytia panela  -> panel_all.csv / ciele na imputaciu
    3) cielena imputacia (targeted_impute.sh, Beagle + 1000G)  [ak nieco chyba]
    4) doplnenie imputovanych genotypov do predispozicneho normalizovaneho suboru

Vysledny predispozicny normalizovany subor sa skopiruje do cieloveho priecinka,
ktory pouzivatel zada v UI (akceptuje aj Windows cestu C:\\... -> /mnt/c/...).

Vzorky bezia SEKVENCNE — targeted_impute.sh pouziva zdielany scratch adresar
work_ti/, takze dva behy naraz by si liezli do cesty. Samotna imputacia je
vnutorne paralelna (JOBS regionov naraz), takze stroj je aj tak vytazeny.
"""

from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import uuid
from pathlib import Path
from typing import List, Optional

from fastapi import APIRouter, File, Form, HTTPException, UploadFile

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dna_tools.parsers import parse_file, relabel_from_reference       # noqa: E402
from dna_tools.convert import write_23andme, merge_imputed_into_23andme  # noqa: E402
from dna_tools.panel import check_panel, annotate_by_rsid              # noqa: E402
from dna_tools.targeted_impute import merge_panel                      # noqa: E402
from dna_tools import wgs                                              # noqa: E402
from app.runtime_security import (                                     # noqa: E402
    JOB_LOCK,
    MAX_BATCH_SAMPLES,
    PUBLIC_MODE,
    prune_runtime_files,
    save_upload_limited,
)

WORK = ROOT / "work"
DATA = ROOT / "data"
LOGS = ROOT / "logs"
ID_REF = DATA / "id_position_reference_v3v4v5.tsv.gz"
DEFAULT_PANEL = DATA / "panel_full_561.csv"
BATCH_LOG = LOGS / "batch.log"

router = APIRouter()

BATCH_JOBS: dict = {}


# --- pomocne ------------------------------------------------------------------

def _blog(msg: str) -> None:
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    try:
        LOGS.mkdir(exist_ok=True)
        with open(BATCH_LOG, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")
    except Exception:  # noqa: BLE001
        pass


_WIN_DRIVE = re.compile(r"^([A-Za-z]):[\\/](.*)$")


def resolve_out_dir(raw: str) -> Path:
    """Prelozi pouzivatelom zadanu cestu na realnu cestu v tomto systeme.

    Prijima:
      - prazdne / relativne  -> ROOT/<raw>  (default 'out')
      - Windows 'C:\\Users\\x' -> '/mnt/c/Users/x'  (ked bezime vo WSL)
      - absolutne POSIX cesty -> tak ako su
    """
    raw = (raw or "").strip().strip('"').strip("'")
    if not raw:
        raw = "out"
    m = _WIN_DRIVE.match(raw)
    if m:
        drive, rest = m.group(1).lower(), m.group(2).replace("\\", "/")
        # vo WSL su Windows disky pod /mnt/c, /mnt/d, ...
        if not Path(f"/mnt/{drive}").is_dir():
            raise ValueError(
                f"Windows cestu '{raw}' sa nedá preložiť — disk {drive.upper()}: "
                f"nie je namountovaný (/mnt/{drive} neexistuje). "
                f"Zadaj radšej relatívnu cestu (napr. 'out') alebo cestu v štýle /mnt/c/…")
        return Path(f"/mnt/{drive}") / rest
    p = Path(raw.replace("\\", "/"))
    if not p.is_absolute():
        p = ROOT / p
    return p


def _resolve_bam_dir(raw: str) -> str:
    """Ako resolve_out_dir, ale nezhodi beh — ked sa neda prelozit, vrati prazdny retazec."""
    raw = (raw or "").strip().strip('"').strip("'")
    if not raw:
        return ""
    try:
        return str(resolve_out_dir(raw))
    except ValueError:
        return raw


def _save_upload(up: UploadFile, prefix: str) -> Path:
    prune_runtime_files(WORK)
    dest = WORK / f"{prefix}_{Path(up.filename).name}"
    return save_upload_limited(up, dest)


def _safe_stem(name: str) -> str:
    stem = Path(name).stem
    return re.sub(r"[^A-Za-z0-9._-]+", "_", stem) or "vzorka"


def find_bam(bam_dir: str, vcf_name: str) -> Optional[str]:
    """Najde BAM patriaci k VCF podla nazvu: S24061.GATK.snp.vcf -> S24061*.bam.

    Skusa postupne kratsie prefixy oddelene bodkou / podtrznikom, aby sedeli
    aj nazvy typu 'vzorka_1.GATK.snp.vcf' <-> 'vzorka_1.sorted.dedup.bam'.
    """
    if not bam_dir:
        return None
    d = Path(bam_dir)
    if not d.is_dir():
        return None
    stem = Path(vcf_name).name
    for suf in (".vcf.gz", ".vcf"):
        if stem.endswith(suf):
            stem = stem[: -len(suf)]
            break
    cands = [stem]
    while "." in cands[-1]:
        cands.append(cands[-1].rsplit(".", 1)[0])
    for base in cands:
        hits = sorted(d.glob(base + "*.bam")) + sorted(d.glob(base + "*.cram"))
        if hits:
            return str(hits[0])
    return None


def _is_vcf(path: Path) -> bool:
    n = str(path).lower()
    return n.endswith(".vcf") or n.endswith(".vcf.gz")


def _imputation_ready() -> bool:
    ref_dir = os.environ.get("REF_DIR", "")
    beagle = os.environ.get("BEAGLE_JAR", "")
    fasta = os.environ.get("FASTA", "")
    bref3 = len(list(Path(ref_dir).glob("*.bref3"))) if ref_dir and Path(ref_dir).is_dir() else 0
    return bool(shutil.which("java") and beagle and Path(beagle).exists()
                and bref3 > 0 and fasta and Path(fasta).exists())


# --- jadro behu jednej vzorky --------------------------------------------------

def _item_set(item: dict, **kw) -> None:
    item.update(kw)


def _wes_plan(src: Path, orig_name: str, job: dict):
    """Ma sa tento vstup spracovat ako WGS/WES? Vrati (bam_path|None, build) alebo None.

    WES vetva sa pouzije ked je vstup VCF a bud je v hg38 (vtedy je liftover
    povinny), alebo k nemu existuje BAM (vtedy vieme doplnit ref/ref). Cipovy
    VCF v GRCh37 bez BAM ide starou cestou — tam sa nic nezlepsi.
    """
    # POZOR: parovat BAM podla PODVODNEHO nazvu uploadu (orig_name), nie podla
    # cesty v work/ — tam ma subor pred nazvom nahodny prefix a glob by nesadol.
    if not _is_vcf(Path(orig_name)):
        return None
    build, why = wgs.detect_vcf_build(str(src))
    bam = find_bam(job["opts"].get("bam_dir", ""), orig_name)
    if build == "hg38" or bam:
        return bam, build
    return None


def _run_impute(job: dict, item: dict, sample_txt: Path, targets: Path,
                status_csv: Path, out_prefix: Path, env: dict) -> Optional[Path]:
    """Spusti targeted_impute.sh a strameuje progres do item. Vrati panel_final.csv."""
    cmd = ["bash", str(ROOT / "targeted_impute.sh"), str(sample_txt), str(targets),
           str(status_csv), str(out_prefix)]
    _blog(f"  impute: {' '.join(cmd)}")
    logf = open(LOGS / "impute.log", "a", encoding="utf-8")
    logf.write(f"\n===== batch {job['id'][:8]} / {item['name']} @ "
               f"{time.strftime('%Y-%m-%d %H:%M:%S')} =====\n")
    logf.flush()
    proc = subprocess.Popen(cmd, cwd=str(ROOT), env=env, text=True, bufsize=1,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    job["proc"] = proc
    for raw in proc.stdout:
        line = raw.rstrip("\n")
        logf.write(line + "\n")
        logf.flush()
        job["tail"].append(f"[{item['name']}] {line}")
        if len(job["tail"]) > 80:
            job["tail"].pop(0)
        m = re.search(r"regionov:\s*(\d+)", line)
        if m:
            item["regions_total"] = int(m.group(1))
        if "REGION_DONE" in line:
            item["regions_done"] = item.get("regions_done", 0) + 1
    proc.wait()
    job["proc"] = None
    logf.write(f"----- exit {proc.returncode} -----\n")
    logf.close()
    if proc.returncode != 0:
        raise RuntimeError(f"targeted_impute.sh skoncil s kodom {proc.returncode} "
                           f"(detaily v logs/impute.log)")
    final = Path(str(out_prefix) + "_panel_final.csv")
    if not final.exists():
        raise RuntimeError("imputacia nevyprodukovala _panel_final.csv")
    return final


def _process_one(job: dict, item: dict, src: Path, panel_path: Path, out_dir: Path) -> None:
    prefix = job["prefix"]
    stem = item["stem"]
    online = job["opts"]["online"]

    # --- 1) prevod -----------------------------------------------------------
    _item_set(item, stage="prevod", stage_no=1, status="running")
    conv = WORK / f"{prefix}_{stem}_normalized_genotype.txt"

    wes = _wes_plan(src, item["name"], job)
    if wes is not None:
        # WGS/WES vetva: liftover hg38->GRCh37 + doplnenie ref/ref podla pokrytia
        bam, build = wes
        _item_set(item, stage=f"prevod WES ({build})", bam=bam)
        _blog(f"  [{stem}] WES/WGS VCF, build={build}, BAM={bam or 'žiadny'}")
        r = wgs.wes_vcf_to_23andme(
            str(src), str(conv), fasta_b37=os.environ.get("FASTA", ""),
            id_reference=str(ID_REF), bam_path=bam,
            min_dp=job["opts"]["min_dp"], sample_name=stem,
            progress=lambda m: (job["tail"].append(f"[{stem}] {m}"),
                                job["tail"].pop(0) if len(job["tail"]) > 80 else None))
        _item_set(item, source_format=f"vcf/{build}", n_variants=r.n_template,
                  n_called=r.n_variant + r.n_homref, converted=conv.name,
                  wes_variant=r.n_variant, wes_homref=r.n_homref, wes_nocall=r.n_nocall,
                  wes_bam=bool(r.bam_used))
        for w in r.warnings:
            job["tail"].append(f"[{stem}] ⚠ {w}")
        _blog(f"  [{stem}] WES prevod OK: {r.n_variant} variantov, {r.n_homref} ref/ref, "
              f"{r.n_nocall} no-call z {r.n_template} šablónových pozícií")
    else:
        data = parse_file(str(src))
        if data.source_format == "illumina_final_report":
            if ID_REF.exists():
                relabel_from_reference(data, str(ID_REF))
            annotate_by_rsid(data, str(ID_REF) if ID_REF.exists() else None, online=online)
        write_23andme(data, str(conv), drop_no_call=False, keep_non_rs=True)
        _item_set(item, source_format=data.source_format, n_variants=data.n_variants,
                  n_called=data.n_called, converted=conv.name)
        _blog(f"  [{stem}] prevod OK ({data.source_format}, {data.n_variants} variantov)")

    # Sibling/IBD výstup je nemenná kópia priamo nameraných genotypov. Nikdy
    # doň neskôr nepridávame Beagle odhady.
    sibling = WORK / f"{prefix}_{stem}_sibling_measured_normalized.txt"
    shutil.copyfile(conv, sibling)

    # --- 2) panel-check ------------------------------------------------------
    _item_set(item, stage="kontrola panela", stage_no=2)
    res = check_panel(str(conv), str(panel_path),
                      id_reference=str(ID_REF) if ID_REF.exists() else None,
                      online=online)
    status_csv = WORK / f"{prefix}_{stem}_panel_all.csv"
    res.table.to_csv(status_csv, index=False)
    targets_csv = WORK / f"{prefix}_{stem}_panel_missing.csv"
    missing = res.table[res.table["status"] != "called"]
    missing.to_csv(targets_csv, index=False)
    _item_set(item, n_panel=res.n_panel, called_before=res.called,
              coverage_before=round(100 * res.coverage, 1), to_impute=int(len(missing)))
    _blog(f"  [{stem}] panel: {res.called}/{res.n_panel} zavolanych "
          f"({100*res.coverage:.1f}%), na imputaciu {len(missing)}")

    # --- 3) imputacia --------------------------------------------------------
    completed = conv
    panel_final = WORK / f"{prefix}_{stem}_predisposition_panel.csv"
    have_targets = any(str(r.get("position", "")).strip() not in ("", "nan", "<NA>", "None")
                       for r in csv.DictReader(open(targets_csv, encoding="utf-8")))
    if len(missing) and have_targets:
        _item_set(item, stage="imputácia (Beagle)", stage_no=3)
        out_prefix = WORK / f"{prefix}_{stem}_targeted"
        env = {**os.environ,
               "DR2MIN": str(job["opts"]["dr2_min"]),
               "DEDUP_BY": job["opts"]["dedup_by"],
               "FLANK": str(job["opts"]["flank"])}
        generated_panel = _run_impute(job, item, conv, targets_csv, status_csv, out_prefix, env)
        shutil.copyfile(generated_panel, panel_final)
        _item_set(item, panel_final=panel_final.name)

        # --- 4) doplnenie do predispozicneho normalizovaneho vystupu --------
        _item_set(item, stage="doplnenie do normalizovaného výstupu", stage_no=4)
        completed = WORK / f"{prefix}_{stem}_predisposition_normalized.txt"
        stats = merge_imputed_into_23andme(str(conv), str(panel_final), str(completed),
                                           min_conf=job["opts"]["min_conf"])
        _item_set(item, replaced_nocall=stats.get("replaced_nocall", 0),
                  added_absent=stats.get("added_absent", 0))
        _blog(f"  [{stem}] merge: +{stats.get('added_absent',0)} doplnenych, "
              f"{stats.get('replaced_nocall',0)} prepisanych no-call")
        # dopocitaj pokrytie po imputacii
        try:
            res2 = check_panel(str(completed), str(panel_path),
                               id_reference=str(ID_REF) if ID_REF.exists() else None,
                               online=False)
            _item_set(item, called_after=res2.called,
                      coverage_after=round(100 * res2.coverage, 1))
        except Exception:  # noqa: BLE001
            pass
    else:
        merge_panel(str(status_csv)).to_csv(panel_final, index=False)
        why = ("panel je kompletne pokryty" if not len(missing)
               else "chybajuce lokusy nemaju znamu poziciu (nedaju sa imputovat)")
        _item_set(item, stage="imputácia preskočená", stage_no=4,
                  skip_reason=("panel je už kompletne pokrytý" if not len(missing)
                               else "chýbajúce lokusy nemajú známu pozíciu"),
                  called_after=item.get("called_before"),
                  coverage_after=item.get("coverage_before"))
        _blog(f"  [{stem}] imputacia preskocena — {why}")

    # --- 5) ulozenie do cieloveho priecinka ---------------------------------
    _item_set(item, stage="ukladám výsledok", stage_no=5)
    out_dir.mkdir(parents=True, exist_ok=True)
    sibling_dest = out_dir / f"{stem}_sibling_measured_normalized.txt"
    predisposition_dest = out_dir / f"{stem}_predisposition_normalized.txt"
    panel_dest = out_dir / f"{stem}_predisposition_panel.csv"
    shutil.copyfile(sibling, sibling_dest)
    shutil.copyfile(completed, predisposition_dest)
    shutil.copyfile(panel_final, panel_dest)
    _item_set(item, stage="hotovo", stage_no=6, status="ok",
              out_path=str(predisposition_dest), out_name=predisposition_dest.name,
              sibling_download=f"/api/download/{sibling.name}",
              predisposition_download=f"/api/download/{Path(completed).name}",
              panel_download=f"/api/download/{panel_final.name}",
              download=f"/api/download/{Path(completed).name}")
    _blog(f"  [{stem}] ULOZENE -> sibling + predisposition + panel")


def _run_batch(job_id: str) -> None:
    job = BATCH_JOBS[job_id]
    job["running"] = True
    job["started"] = time.time()
    out_dir = Path(job["out_dir_resolved"])
    panel_path = Path(job["panel_path"])
    _blog(f"=== DAVKA {job_id[:8]}: {len(job['items'])} vzoriek -> {out_dir} ===")
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:  # noqa: BLE001
        job["error"] = f"Cieľový priečinok sa nedá vytvoriť: {out_dir} ({e})"
        job["running"] = False
        job["done"] = True
        job["ended"] = time.time()
        return

    for item in job["items"]:
        if job.get("cancel"):
            _item_set(item, stage="zrušené", status="cancelled")
            continue
        item["started"] = time.time()
        try:
            _process_one(job, item, Path(item["src"]), panel_path, out_dir)
        except Exception as e:  # noqa: BLE001
            _item_set(item, stage="chyba", status="error", error=str(e))
            _blog(f"  CHYBA [{item['name']}]: {e}")
        finally:
            item["ended"] = time.time()
            job["done_count"] = sum(1 for it in job["items"]
                                    if it.get("status") in ("ok", "error", "cancelled"))

    # --- suhrn dávky ---------------------------------------------------------
    try:
        rep = out_dir / "_prehlad_davky.csv"
        cols = ["vzorka", "format", "variantov", "panel_lokusov",
                "pokrytie_pred_%", "pokrytie_po_%", "doimputovanych",
                "vysledny_subor", "stav", "chyba"]
        with open(rep, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(cols)
            for it in job["items"]:
                w.writerow([
                    it.get("name", ""), it.get("source_format", ""),
                    it.get("n_variants", ""), it.get("n_panel", ""),
                    it.get("coverage_before", ""), it.get("coverage_after", ""),
                    it.get("added_absent", ""), it.get("out_name", ""),
                    it.get("status", ""), it.get("error", ""),
                ])
        job["report"] = str(rep)
    except Exception as e:  # noqa: BLE001
        _blog(f"  prehlad davky sa nepodarilo zapisat: {e}")

    job["running"] = False
    job["done"] = True
    job["ended"] = time.time()
    ok = sum(1 for it in job["items"] if it.get("status") == "ok")
    _blog(f"=== DAVKA {job_id[:8]} HOTOVA: {ok}/{len(job['items'])} OK, "
          f"{int(job['ended']-job['started'])} s ===")


def _run_batch_locked(job_id: str) -> None:
    try:
        _run_batch(job_id)
    finally:
        if JOB_LOCK.locked():
            JOB_LOCK.release()


# --- endpointy ----------------------------------------------------------------

@router.post("/api/batch-run")
async def api_batch_run(samples: List[UploadFile] = File(...),
                        panel: Optional[UploadFile] = File(None),
                        out_dir: str = Form("out"),
                        online: bool = Form(False),
                        dr2_min: float = Form(0.9),
                        flank: int = Form(250000),
                        dedup_by: str = Form("position"),
                        min_conf: float = Form(0.9),
                        bam_dir: str = Form(""),
                        min_dp: int = Form(8)):
    """Spusti celu pipeline pre vsetky nahrane vzorky. Vrati job_id."""
    if not samples:
        raise HTTPException(400, "Nevybral si žiadnu vzorku.")
    if len(samples) > MAX_BATCH_SAMPLES:
        raise HTTPException(400, f"Naraz možno spracovať najviac {MAX_BATCH_SAMPLES} vzoriek.")
    if not _imputation_ready():
        raise HTTPException(400,
            "Imputácia nie je pripravená (chýba java / BEAGLE_JAR / REF_DIR / FASTA). "
            "Spusti appku cez '▶ Genome Converter' alebo `bash start_app.sh` vo WSL.")

    if not JOB_LOCK.acquire(blocking=False):
        raise HTTPException(409, "Iná imputácia už beží. Počkaj na jej dokončenie.")

    prefix = uuid.uuid4().hex[:8]
    try:
        if PUBLIC_MODE:
            resolved = WORK / f"{prefix}_results"
            bam_dir = ""
            online = False
        else:
            resolved = resolve_out_dir(out_dir)
        if panel is None:
            if not DEFAULT_PANEL.exists():
                raise HTTPException(500, "Vstavaný 561-SNP panel nie je dostupný.")
            panel_path = DEFAULT_PANEL
            panel_name = DEFAULT_PANEL.name
        else:
            panel_path = _save_upload(panel, prefix)
            panel_name = panel.filename

        items = []
        for up in samples:
            src = _save_upload(up, prefix)
            items.append({
                "name": up.filename, "stem": _safe_stem(up.filename), "src": str(src),
                "stage": "čaká", "stage_no": 0, "status": "pending", "error": None,
                "regions_done": 0, "regions_total": 0,
                "started": None, "ended": None, "out_path": None, "out_name": None,
            })
    except ValueError as e:
        JOB_LOCK.release()
        raise HTTPException(400, str(e))
    except Exception:
        JOB_LOCK.release()
        raise

    job_id = uuid.uuid4().hex
    BATCH_JOBS[job_id] = {
        "id": job_id, "prefix": prefix, "items": items,
        "panel_path": str(panel_path), "panel_name": panel_name,
        "out_dir_input": out_dir, "out_dir_resolved": str(resolved),
        "opts": {"online": bool(online), "dr2_min": float(dr2_min), "flank": int(flank),
                 "dedup_by": "rsid" if str(dedup_by).lower().startswith("rs") else "position",
                 "min_conf": float(min_conf),
                 "bam_dir": _resolve_bam_dir(bam_dir),
                 "min_dp": int(min_dp)},
        "tail": [], "running": False, "done": False, "done_count": 0,
        "error": None, "report": None, "started": None, "ended": None,
        "cancel": False, "proc": None,
    }
    try:
        threading.Thread(target=_run_batch_locked, args=(job_id,), daemon=True).start()
    except Exception:
        BATCH_JOBS.pop(job_id, None)
        JOB_LOCK.release()
        raise
    return {"job_id": job_id, "out_dir": str(resolved), "n": len(items)}


@router.get("/api/batch-progress/{job_id}")
def api_batch_progress(job_id: str):
    job = BATCH_JOBS.get(job_id)
    if job is None:
        raise HTTPException(404, "neznáma dávka")
    now = time.time()
    elapsed = int((job["ended"] or now) - job["started"]) if job["started"] else 0
    items = []
    for it in job["items"]:
        el = 0
        if it.get("started"):
            el = int((it.get("ended") or now) - it["started"])
        items.append({
            "name": it["name"], "stage": it["stage"], "stage_no": it["stage_no"],
            "status": it["status"], "error": it["error"],
            "regions_done": it.get("regions_done", 0),
            "regions_total": it.get("regions_total", 0),
            "source_format": it.get("source_format"), "n_variants": it.get("n_variants"),
            "n_panel": it.get("n_panel"),
            "coverage_before": it.get("coverage_before"),
            "coverage_after": it.get("coverage_after"),
            "added_absent": it.get("added_absent"),
            "skip_reason": it.get("skip_reason"),
            "wes_variant": it.get("wes_variant"), "wes_homref": it.get("wes_homref"),
            "wes_nocall": it.get("wes_nocall"), "wes_bam": it.get("wes_bam"),
            "bam": it.get("bam"),
            "out_name": it.get("out_name"), "out_path": it.get("out_path"),
            "download": it.get("download"),
            "sibling_download": it.get("sibling_download"),
            "predisposition_download": it.get("predisposition_download"),
            "panel_download": it.get("panel_download"),
            "elapsed": el,
        })
    return {
        "out_dir": job["out_dir_resolved"], "panel": job["panel_name"],
        "running": job["running"], "done": job["done"], "error": job["error"],
        "done_count": job["done_count"], "total": len(job["items"]),
        "report": job["report"], "elapsed": elapsed,
        "tail": job["tail"][-20:], "items": items,
    }


@router.post("/api/batch-cancel/{job_id}")
def api_batch_cancel(job_id: str):
    job = BATCH_JOBS.get(job_id)
    if job is None:
        raise HTTPException(404, "neznáma dávka")
    job["cancel"] = True
    proc = job.get("proc")
    if proc is not None:
        try:
            proc.terminate()
        except Exception:  # noqa: BLE001
            pass
    _blog(f"davka {job_id[:8]}: ZRUSENA pouzivatelom")
    return {"ok": True}
