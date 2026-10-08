"""FastAPI UI for local or password-protected hosted genome normalization.

Spustenie (lokalne):
    cd genome-app
    pip install -r requirements.txt
    uvicorn app.main:app --host 127.0.0.1 --port 8000
    -> otvor http://127.0.0.1:8000

Pri lokálnom spustení dáta neopúšťajú počítač. Verejný režim musí mať
nastavené PUBLIC_MODE=1, APP_USERNAME a APP_PASSWORD.
"""

from __future__ import annotations

import csv
import io
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import uuid
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, File, Form, UploadFile, HTTPException, Request
from fastapi.responses import HTMLResponse, FileResponse, PlainTextResponse
from fastapi.staticfiles import StaticFiles

# --- cesty --------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent   # genome-app/
sys.path.insert(0, str(ROOT))                    # aby sa importol balik dna_tools

from dna_tools.detect import detect_format_verbose            # noqa: E402
from dna_tools.parsers import parse_file, relabel_from_reference   # noqa: E402
from dna_tools.convert import write_23andme, merge_imputed_into_23andme   # noqa: E402
from dna_tools.panel import check_panel, annotate_by_rsid     # noqa: E402
from dna_tools.targeted_impute import build_windows           # noqa: E402
from dna_tools.analysis_prep import prepare_for_analyses      # noqa: E402
from dna_tools import wgs                                     # noqa: E402
from app.runtime_security import (                            # noqa: E402
    JOB_LOCK,
    MAX_REQUEST_MB,
    PUBLIC_MODE,
    RETENTION_HOURS,
    basic_auth_role,
    prune_runtime_files,
    save_upload_limited,
    test_quota_status,
    validate_public_config,
)

APP_VERSION = "2.4"

WORK = ROOT / "work"
LOGS = ROOT / "logs"
DATA = ROOT / "data"
STATIC = ROOT / "static"
DEMO = ROOT / "demo"
for d in (WORK, LOGS, DATA):
    d.mkdir(exist_ok=True)

LOG_FILE = LOGS / "app.log"
IMPUTE_LOG = LOGS / "impute.log"   # zivy vypis behu imputacie (tail -f)
DEMO_LOCK = threading.Lock()

DEMO_ALLOWED_PATHS = {
    "/",
    "/api/status",
    "/api/demo-run",
}


def demo_path_allowed(path: str) -> bool:
    """Return whether the shared DEMO role may access an HTTP path."""

    return (
        path in DEMO_ALLOWED_PATHS
        or path.startswith("/api/demo-download/")
        or path.startswith("/static/")
    )


def test_path_allowed(path: str) -> bool:
    """Return whether the shared test role may use an HTTP path."""

    return (
        path in {"/", "/api/status", "/api/batch-run", "/api/demo-run"}
        or path.startswith("/api/batch-progress/")
        or path.startswith("/api/batch-cancel/")
        or path.startswith("/api/download/")
        or path.startswith("/api/demo-download/")
        or path.startswith("/static/")
    )


def log(msg: str) -> None:
    """Zapise riadok do logu (aj na stdout). Log si moze citat kto ladi appku."""
    line = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a", encoding="utf-8") as fh:
        fh.write(line + "\n")


validate_public_config()
app = FastAPI(
    title="Genome Normalizer",
    version=APP_VERSION,
    docs_url=None if PUBLIC_MODE else "/docs",
    redoc_url=None if PUBLIC_MODE else "/redoc",
    openapi_url=None if PUBLIC_MODE else "/openapi.json",
)
log(f"=== Genome Normalizer v{APP_VERSION} startuje ===")


@app.middleware("http")
async def hosted_security(request, call_next):
    """Require credentials in public mode and add conservative browser headers."""

    role = "owner"
    if PUBLIC_MODE and request.url.path != "/healthz":
        role = basic_auth_role(request.headers.get("authorization"))
        if role is None:
            return PlainTextResponse(
                "Vyžaduje sa prihlásenie.",
                status_code=401,
                headers={"WWW-Authenticate": 'Basic realm="Genome Normalizer"'},
            )
        if role == "demo" and not demo_path_allowed(request.url.path):
            return PlainTextResponse(
                "DEMO účet povoľuje iba zabudovanú syntetickú ukážku.",
                status_code=403,
            )
        if role == "test" and not test_path_allowed(request.url.path):
            return PlainTextResponse(
                "Testovací účet povoľuje kompletnú dávkovú pipeline do denného limitu.",
                status_code=403,
            )
    request.state.account_role = role
    content_length = request.headers.get("content-length")
    if request.method in {"POST", "PUT", "PATCH"} and content_length:
        try:
            if int(content_length) > MAX_REQUEST_MB * 1024 * 1024:
                return PlainTextResponse("Požiadavka je príliš veľká.", status_code=413)
        except ValueError:
            pass
    response = await call_next(request)
    response.headers["Cache-Control"] = "no-store"
    response.headers["X-Content-Type-Options"] = "nosniff"
    response.headers["X-Frame-Options"] = "DENY"
    response.headers["Referrer-Policy"] = "no-referrer"
    response.headers["Permissions-Policy"] = "camera=(), microphone=(), geolocation=()"
    if PUBLIC_MODE:
        response.headers["Strict-Transport-Security"] = "max-age=31536000; includeSubDomains"
    return response

# davkove spracovanie (jedno tlacidlo -> cela pipeline pre N vzoriek)
from app.batch import router as batch_router   # noqa: E402
app.include_router(batch_router)


def _cleanup_loop() -> None:
    interval = max(300, min(1800, int(RETENTION_HOURS * 1800)))
    while True:
        time.sleep(interval)
        if not JOB_LOCK.locked():
            prune_runtime_files(WORK)


@app.on_event("startup")
def start_runtime_cleanup() -> None:
    if PUBLIC_MODE:
        prune_runtime_files(WORK)
        threading.Thread(target=_cleanup_loop, name="runtime-cleanup", daemon=True).start()


# --- pomocne ------------------------------------------------------------------

def _save_upload(up: UploadFile) -> Path:
    prune_runtime_files(WORK)
    sid = uuid.uuid4().hex[:8]
    dest = WORK / f"{sid}_{Path(up.filename).name}"
    save_upload_limited(up, dest)
    display_name = dest.name if PUBLIC_MODE else up.filename
    log(f"upload: {display_name} ({dest.stat().st_size} B)")
    return dest


def _which(name: str) -> Optional[str]:
    return shutil.which(name)


# --- endpointy ----------------------------------------------------------------

@app.get("/healthz")
def healthz():
    return {"ok": True, "version": APP_VERSION}

@app.get("/", response_class=HTMLResponse)
def index():
    idx = STATIC / "index.html"
    if idx.exists():
        return idx.read_text(encoding="utf-8")
    return "<h1>Genome Normalizer</h1><p>static/index.html chýba.</p>"


@app.get("/api/status")
def status(request: Request):
    """Zisti dostupnost nastrojov a referencie (na imputaciu)."""
    ref_dir = os.environ.get("REF_DIR", "")
    beagle = os.environ.get("BEAGLE_JAR", "")
    fasta = os.environ.get("FASTA", "")
    bref3 = 0
    if ref_dir and Path(ref_dir).is_dir():
        bref3 = len(list(Path(ref_dir).glob("*.bref3")))
    payload = {
        "version": APP_VERSION,
        "public_mode": PUBLIC_MODE,
        "account_role": getattr(request.state, "account_role", "owner"),
        "demo_available": (DEMO / "demo_myheritage_raw.csv").exists(),
        "retention_hours": RETENTION_HOURS,
        "default_panel": (DATA / "panel_full_561.csv").exists(),
        "python": sys.version.split()[0],
        "java": bool(_which("java")),
        "beagle_jar": bool(beagle and Path(beagle).exists()),
        "fasta": bool(fasta and Path(fasta).exists()),
        "bref3_files": bref3,
        "id_reference": (DATA / "id_position_reference_v3v4v5.tsv.gz").exists(),
        "imputation_ready": bool(_which("java") and beagle and Path(beagle or "x").exists()
                                 and bref3 > 0 and fasta and Path(fasta or "x").exists()),
        "wgs": wgs.wgs_status(),
    }
    if getattr(request.state, "account_role", "owner") == "test":
        payload["test_quota"] = test_quota_status()
    return payload


@app.post("/api/detect")
async def api_detect(file: UploadFile = File(...)):
    path = _save_upload(file)
    try:
        fmt, reason = detect_format_verbose(str(path))
        log(f"detect: {path.name} -> {fmt.value} ({reason})")
        return {"format": fmt.value, "reason": reason, "token": path.name}
    except Exception as e:  # noqa: BLE001
        log(f"ERROR detect: {e}")
        raise HTTPException(500, str(e))


@app.post("/api/convert")
async def api_convert(file: UploadFile = File(...),
                      drop_no_call: bool = Form(False),
                      only_rs: bool = Form(False),
                      online: bool = Form(False)):
    path = _save_upload(file)
    try:
        data = parse_file(str(path))
        # Illumina: prevezmi kanonicke oznacenie podla lokalnej pozicnej referencie.
        if data.source_format == "illumina_final_report":
            id_ref = DATA / "id_position_reference_v3v4v5.tsv.gz"
            if id_ref.exists():
                n = relabel_from_reference(data, str(id_ref))
                data.warnings.append(
                    "Illumina: %d markerov preznacenych podla lokalnej ID referencie (podla pozicie); "
                    "zvysok ponechany s povodnym oznacenim." % n)
            # rs-only sondy (chr 0, ale maju rsID) — dopln suradnice podla rs cisla:
            # lokalne z ID referencie, a ak online=True, zvysok cez Ensembl GRCh37
            na = annotate_by_rsid(data, str(id_ref) if id_ref.exists() else None, online=online)
            if na:
                data.warnings.append(
                    "Illumina: %d rs-only sondam doplnene chrom/pozicia podla rsID "
                    "(%s)." % (na, "lokal+Ensembl" if online else "lokalna referencia"))
        if data.source_format == "vcf":
            build, why = wgs.detect_vcf_build(str(path))
            if build == "hg38":
                data.warnings.insert(0,
                    "POZOR: tento VCF je v hg38 (%s), ale appka pracuje v GRCh37 — pozicie "
                    "sa nebudu zhodovat a panel vyjde 0%%. Pouzi kartu 'Cela pipeline jednym "
                    "tlacidlom' hore, ktora hg38 VCF prelozi na GRCh37 (a s BAM doplni ref/ref)." % why)
            elif build == "unknown":
                data.warnings.insert(0,
                    "Build tohto VCF sa neda urcit z hlavicky — over, ze su suradnice v GRCh37.")
        out = WORK / (path.stem + "_normalized_genotype.txt")
        write_23andme(data, str(out), drop_no_call=drop_no_call, keep_non_rs=not only_rs)
        log(f"convert: {path.name} ({data.source_format}) -> {out.name} "
            f"[{data.n_variants} variantov, {data.n_called} zavolanych]")
        # nahlad prvych par datovych riadkov vystupu (bez komentarovej hlavicky)
        preview = []
        with open(out, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                preview.append(line.rstrip("\n"))
                if len(preview) >= 20:
                    break
        return {
            "source_format": data.source_format,
            "n_variants": data.n_variants,
            "n_called": data.n_called,
            "no_call": data.n_variants - data.n_called,
            "warnings": data.warnings,
            "preview": preview,
            "download": f"/api/download/{out.name}",
            "output": out.name,
        }
    except Exception as e:  # noqa: BLE001
        log(f"ERROR convert: {e}")
        raise HTTPException(500, str(e))


@app.post("/api/panel-check")
async def api_panel_check(sample: UploadFile = File(...),
                          panel: UploadFile = File(...),
                          online: bool = Form(False)):
    sp = _save_upload(sample)
    pp = _save_upload(panel)
    try:
        id_ref = DATA / "id_position_reference_v3v4v5.tsv.gz"
        res = check_panel(str(sp), str(pp),
                          id_reference=str(id_ref) if id_ref.exists() else None,
                          online=online)
        out = WORK / (sp.stem + "_panel_all.csv")
        res.table.to_csv(out, index=False)
        miss = WORK / (sp.stem + "_panel_missing.csv")
        res.table[res.table["status"] != "called"].to_csv(miss, index=False)
        log(f"panel-check: {sp.name} vs {pp.name} -> called {res.called}/{res.n_panel}"
            f" (doplnenych: ref {res.enriched}, online {res.enriched_online})")
        return {
            "n_panel": res.n_panel, "called": res.called, "no_call": res.no_call,
            "absent": res.absent, "coverage_pct": round(100 * res.coverage, 1),
            "matched_by_rsid": res.matched_by_rsid, "matched_by_pos": res.matched_by_pos,
            "enriched": res.enriched, "enriched_online": res.enriched_online,
            "download_all": f"/api/download/{out.name}",
            "download_missing": f"/api/download/{miss.name}",
        }
    except Exception as e:  # noqa: BLE001
        log(f"ERROR panel-check: {e}")
        raise HTTPException(500, str(e))


@app.post("/api/impute-plan")
async def api_impute_plan(targets: UploadFile = File(...), flank: int = Form(250000)):
    """Offline cast: vypocet okien okolo cielov (nevyzaduje panel)."""
    tp = _save_upload(targets)
    try:
        rows = [(r["chromosome"], int(float(r["position"])))
                for r in csv.DictReader(open(tp)) if r.get("position")]
        regs = build_windows(rows, flank=flank)
        total = sum(r.end - r.start for r in regs)
        log(f"impute-plan: {len(rows)} cielov -> {len(regs)} regionov, {total/1e6:.1f} Mb")
        return {"targets": len(rows), "regions": len(regs),
                "covered_mb": round(total / 1e6, 1)}
    except Exception as e:  # noqa: BLE001
        log(f"ERROR impute-plan: {e}")
        raise HTTPException(500, str(e))


# --- imputacia bezi na POZADI; UI si tiaha progres cez /api/impute-progress ---
IMPUTE_JOBS: dict = {}          # job_id -> stav behu
_STEP_NAMES = {
    1: "vzorka → VCF", 2: "sort + index", 3: "výpočet okien",
    4: "imputácia (Beagle)", 5: "spájanie regiónov", 6: "extrakcia + panel",
}


def _run_impute_job(job_id: str, cmd, env, final_path: str) -> None:
    """Spusti pipeline a strameuje jeho vystup do IMPUTE_JOBS[job_id]."""
    job = IMPUTE_JOBS[job_id]
    job["running"] = True
    job["started"] = time.time()
    try:
        proc = subprocess.Popen(cmd, cwd=str(ROOT), env=env, text=True, bufsize=1,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        logf = open(IMPUTE_LOG, "a", encoding="utf-8")
        logf.write(f"\n===== impute-run {job_id[:8]} @ {time.strftime('%Y-%m-%d %H:%M:%S')} =====\n")
        logf.flush()
        for raw in proc.stdout:                      # merged stdout+stderr
            line = raw.rstrip("\n")
            logf.write(line + "\n"); logf.flush()    # zivy log pre `tail -f logs/impute.log`
            job["tail"].append(line)
            if len(job["tail"]) > 60:
                job["tail"].pop(0)
            m = re.match(r"\[(\d)/6\]", line)
            if m:
                job["step"] = int(m.group(1))
                job["step_name"] = _STEP_NAMES.get(job["step"], "")
            mt = re.search(r"regionov:\s*(\d+)", line)
            if mt:
                job["regions_total"] = int(mt.group(1))
            # jedna dokoncena oblast (uspech/resume/skip/fail) — paralelny beh
            if "REGION_DONE" in line or line.strip().endswith("jar finished"):
                job["regions_done"] = job.get("regions_done", 0) + 1
        proc.wait()
        try:
            logf.write(f"----- exit {proc.returncode} -----\n"); logf.close()
        except Exception:  # noqa: BLE001
            pass
        job["returncode"] = proc.returncode
        if proc.returncode == 0 and Path(final_path).exists():
            job["download"] = f"/api/download/{Path(final_path).name}"
            unimp = Path(str(final_path).replace("_panel_final.csv", "_unimputed.csv"))
            if unimp.exists():
                job["download_unimputed"] = f"/api/download/{unimp.name}"
        log(f"impute-run[{job_id[:8]}]: exit {proc.returncode}")
    except Exception as e:  # noqa: BLE001
        job["error"] = str(e)
        log(f"ERROR impute-run[{job_id[:8]}]: {e}")
    finally:
        job["running"] = False
        job["done"] = True
        job["ended"] = time.time()
        if JOB_LOCK.locked():
            JOB_LOCK.release()


@app.post("/api/impute-run")
async def api_impute_run(sample: UploadFile = File(...),
                         targets: UploadFile = File(...),
                         panel_status: UploadFile = File(...),
                         dr2_min: float = Form(0.9),
                         dedup_by: str = Form("position")):
    """Spusti cielenu imputaciu na pozadi a vrati job_id (progres cez /api/impute-progress)."""
    st = {
        "imputation_ready": bool(
            _which("java")
            and os.environ.get("BEAGLE_JAR")
            and Path(os.environ.get("BEAGLE_JAR", "x")).exists()
            and os.environ.get("REF_DIR")
            and Path(os.environ.get("REF_DIR", "x")).is_dir()
            and os.environ.get("FASTA")
            and Path(os.environ.get("FASTA", "x")).exists()
        )
    }
    if not st["imputation_ready"]:
        raise HTTPException(400, "Imputacia nie je pripravena: chyba java/beagle/panel/FASTA. "
                                 "Nastav BEAGLE_JAR, REF_DIR, FASTA a stiahni panel (setup_reference.sh).")
    if not JOB_LOCK.acquire(blocking=False):
        raise HTTPException(409, "Iná imputácia už beží. Počkaj na jej dokončenie.")
    try:
        sp = _save_upload(sample)
        tp = _save_upload(targets)
        pp = _save_upload(panel_status)
    except Exception:
        JOB_LOCK.release()
        raise
    out_prefix = WORK / (sp.stem + "_targeted")
    final = Path(str(out_prefix) + "_panel_final.csv")
    cmd = ["bash", str(ROOT / "targeted_impute.sh"), str(sp), str(tp), str(pp), str(out_prefix)]
    job_id = uuid.uuid4().hex
    IMPUTE_JOBS[job_id] = {
        "step": 0, "step_name": "spúšťam…", "regions_done": 0, "regions_total": 0,
        "tail": [], "running": False, "done": False, "returncode": None,
        "download": None, "download_unimputed": None, "error": None,
        "started": None, "ended": None,
    }
    log(f"impute-run[{job_id[:8]}]: spustam {' '.join(cmd)}")
    dedup = "rsid" if str(dedup_by).lower().startswith("rs") else "position"
    try:
        threading.Thread(
            target=_run_impute_job,
            args=(job_id, cmd,
                  {**os.environ, "DR2MIN": str(dr2_min), "DEDUP_BY": dedup}, str(final)),
            daemon=True,
        ).start()
    except Exception:
        JOB_LOCK.release()
        raise
    return {"job_id": job_id}


@app.get("/api/impute-progress/{job_id}")
def api_impute_progress(job_id: str):
    """Vrati aktualny stav behu imputacie."""
    job = IMPUTE_JOBS.get(job_id)
    if job is None:
        raise HTTPException(404, "neznamy job")
    elapsed = 0
    if job["started"]:
        elapsed = int((job["ended"] or time.time()) - job["started"])
    return {
        "step": job["step"], "step_name": job["step_name"],
        "regions_done": job["regions_done"], "regions_total": job["regions_total"],
        "running": job["running"], "done": job["done"], "returncode": job["returncode"],
        "download": job["download"], "download_unimputed": job["download_unimputed"],
        "error": job["error"], "elapsed": elapsed,
        "tail": job["tail"][-18:],
    }


@app.post("/api/merge-final")
async def api_merge_final(sample: UploadFile = File(...),
                          panel_final: UploadFile = File(...),
                          min_conf: float = Form(0.0)):
    """Krok 4: doplní imputované genotypy do normalizovaného súboru z kroku 1."""
    sp = _save_upload(sample)
    pf = _save_upload(panel_final)
    out = WORK / (sp.stem + "_predisposition_normalized.txt")
    try:
        stats = merge_imputed_into_23andme(str(sp), str(pf), str(out), min_conf=min_conf)
        log(f"merge-final: {sp.name} + {pf.name} -> {out.name} "
            f"[prepisanych no-call {stats['replaced_nocall']}, pridanych {stats['added_absent']}]")
        return {**stats, "download": f"/api/download/{out.name}", "output": out.name}
    except Exception as e:  # noqa: BLE001
        log(f"ERROR merge-final: {e}")
        raise HTTPException(500, str(e))


def _demo_downloads() -> dict[str, tuple[Path, str]]:
    demo_work = WORK / "demo"
    return {
        "input": (DEMO / "demo_myheritage_raw.csv", "demo_myheritage_raw.csv"),
        "normalized": (demo_work / "demo_normalized_genotype.txt", "demo_normalized_genotype.txt"),
        "sibling": (demo_work / "demo_sibling_measured_normalized.txt", "demo_sibling_measured_normalized.txt"),
        "panel": (demo_work / "demo_predisposition_panel.csv", "demo_predisposition_panel.csv"),
        "report": (demo_work / "demo_analysis_report.json", "demo_analysis_report.json"),
    }


@app.post("/api/demo-run")
def api_demo_run():
    """Run the real lightweight normalization/QC path on a synthetic sample."""

    source = DEMO / "demo_myheritage_raw.csv"
    if not source.is_file():
        raise HTTPException(503, "Demo vzorka nie je v nasadení dostupná.")
    with DEMO_LOCK:
        demo_work = WORK / "demo"
        demo_work.mkdir(parents=True, exist_ok=True)
        try:
            result = prepare_for_analyses(
                str(source),
                str(demo_work),
                panel_path=str(DATA / "panel_full_561.csv"),
                id_reference_path=str(DATA / "id_position_reference_v3v4v5.tsv.gz"),
                online_annotation=False,
                run_imputation=False,
            )
            normalized = demo_work / "demo_normalized_genotype.txt"
            sibling = demo_work / "demo_sibling_measured_normalized.txt"
            panel = demo_work / "demo_predisposition_panel.csv"
            shutil.copyfile(result.normalized_dataset, normalized)
            shutil.copyfile(result.sibling_dataset, sibling)
            shutil.copyfile(result.predisposition_panel, panel)

            preview = []
            with normalized.open("r", encoding="utf-8") as handle:
                for line in handle:
                    if line.startswith("#") or not line.strip():
                        continue
                    preview.append(line.rstrip("\n").split("\t"))
                    if len(preview) == 10:
                        break

            demo_report = {
                "demo": True,
                "synthetic_data": True,
                "source_file": source.name,
                "source_format": result.source_format,
                "reference_build": "GRCh37",
                "normalized_schema": "rsid, chromosome, position, genotype; plus strand",
                "measured_variants": result.measured_variants,
                "measured_called": result.measured_called,
                "sibling_markers": result.sibling_markers,
                "sibling_contains_imputed_genotypes": False,
                "panel_markers": result.panel_markers,
                "panel_measured": result.panel_measured,
                "panel_imputed": 0,
                "panel_missing_after": result.panel_missing_after,
                "imputation_run": False,
                "note": (
                    "DEMO používa syntetickú vzorku a reálnu normalizáciu s panelovou kontrolou. "
                    "Nákladná Beagle imputácia je pre zdieľané demo konto zámerne vypnutá."
                ),
            }
            report_path = demo_work / "demo_analysis_report.json"
            report_path.write_text(
                json.dumps(demo_report, ensure_ascii=False, indent=2), encoding="utf-8"
            )
            log("demo-run: syntetická MyHeritage vzorka spracovaná")
            return {
                **demo_report,
                "coverage_pct": round(
                    100 * result.panel_measured / result.panel_markers, 1
                ) if result.panel_markers else 0.0,
                "preview": preview,
                "downloads": {
                    key: f"/api/demo-download/{key}" for key in _demo_downloads()
                },
            }
        except Exception as exc:  # noqa: BLE001
            log(f"ERROR demo-run: {exc}")
            raise HTTPException(500, f"Demo analýza zlyhala: {exc}") from exc


@app.get("/api/demo-download/{kind}")
def demo_download(kind: str):
    entry = _demo_downloads().get(kind)
    if entry is None:
        raise HTTPException(404, "Neznámy demo výstup.")
    path, filename = entry
    if not path.is_file():
        raise HTTPException(404, "Najprv spustite demo analýzu.")
    return FileResponse(str(path), filename=filename)


@app.get("/api/download/{name}")
def download(name: str):
    p = (WORK / Path(name).name).resolve()
    if p.parent != WORK.resolve() or not p.is_file():
        raise HTTPException(404, "subor neexistuje")
    return FileResponse(str(p), filename=name)


@app.get("/api/logs", response_class=PlainTextResponse)
def get_logs(tail: int = 200):
    if PUBLIC_MODE:
        raise HTTPException(404, "log nie je vo verejnom režime dostupný")
    if not LOG_FILE.exists():
        return "(log je prazdny)"
    lines = LOG_FILE.read_text(encoding="utf-8").splitlines()
    return "\n".join(lines[-tail:])


# statika (ak by sme chceli servovat aj ine subory)
if STATIC.exists():
    app.mount("/static", StaticFiles(directory=str(STATIC)), name="static")
