"""Prepare one chip sample for both kinship and predisposition analyses.

The pipeline intentionally creates two derived products from one normalized
sample:

* the sibling/IBD dataset contains measured genotypes only;
* the predisposition dataset may contain targeted Beagle imputations.

Keeping the products separate prevents a small, clinically oriented imputation
panel from being mistaken for measured evidence by the IBD analysis.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

import pandas as pd

from .convert import merge_imputed_into_23andme, write_23andme
from .common import NO_CALL
from .impute import build_scaffold
from .panel import annotate_by_rsid, check_panel
from .parsers import parse_file, relabel_from_reference
from .targeted_impute import merge_panel


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_PANEL = ROOT / "data" / "panel_full_561.csv"
DEFAULT_ID_REFERENCE = ROOT / "data" / "id_position_reference_v3v4v5.tsv.gz"


@dataclass
class PreparationResult:
    """Paths and headline QC values produced for one sample."""

    sample: str
    source_format: str
    canonical_profile: str
    normalized_dataset: str
    sibling_dataset: str
    predisposition_panel: str
    predisposition_dataset: Optional[str]
    panel_status: str
    panel_missing: str
    report: str
    measured_variants: int
    measured_called: int
    sibling_markers: int
    panel_markers: int
    panel_measured: int
    panel_imputed: int
    panel_rejected_low_conf: int
    panel_missing_after: int
    imputation_run: bool


def _safe_stem(path: str) -> str:
    stem = Path(path).name
    while Path(stem).suffix.lower() in {".gz", ".zip", ".txt", ".csv", ".tsv", ".vcf"}:
        stem = Path(stem).stem
    return re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("._") or "sample"


def _resolve(path: Optional[str], default: Path) -> Path:
    return Path(path).expanduser().resolve() if path else default.resolve()


def _panel_counts(path: Path) -> dict:
    table = pd.read_csv(path, dtype=str).fillna("")
    source = table["source"].value_counts().to_dict() if "source" in table else {}
    return {
        "measured": int(source.get("measured", 0)),
        "imputed": int(source.get("imputed", 0)),
        "missing": int(source.get("missing", 0)),
    }


def _filter_imputed_panel(source: Path, destination: Path, min_gp: float) -> int:
    """Keep only imputed calls whose per-sample max(GP) reaches the threshold."""

    table = pd.read_csv(source, dtype=str, keep_default_na=False)
    if "source" not in table or "confidence" not in table:
        shutil.copyfile(source, destination)
        return 0
    confidence = pd.to_numeric(table["confidence"], errors="coerce")
    rejected = (table["source"] == "imputed") & (
        confidence.isna() | (confidence < float(min_gp))
    )
    table.loc[rejected, "source"] = "missing"
    table.loc[rejected, "genotype"] = NO_CALL
    table.loc[rejected, "confidence"] = ""
    table.to_csv(destination, index=False)
    return int(rejected.sum())


def _run_targeted_imputation(
    normalized: Path,
    targets: Path,
    status: Path,
    out_prefix: Path,
    work_dir: Path,
    min_gp: float,
    flank: int,
    beagle_jar: Optional[str],
    ref_dir: Optional[str],
    map_dir: Optional[str],
    fasta: Optional[str],
) -> Path:
    """Run the existing Beagle pipeline in a sample-specific scratch folder."""

    script = ROOT / "targeted_impute.sh"
    bash = shutil.which("bash")
    java = shutil.which("java")
    bcftools = shutil.which("bcftools")
    bgzip = shutil.which("bgzip")
    missing_tools = [name for name, value in (
        ("bash", bash), ("java", java), ("bcftools", bcftools), ("bgzip", bgzip)
    ) if not value]
    if missing_tools:
        raise RuntimeError("Chýbajú externé nástroje pre imputáciu: " + ", ".join(missing_tools))

    resolved = {
        "BEAGLE_JAR": _resolve(beagle_jar or os.environ.get("BEAGLE_JAR"), ROOT / "ref" / "beagle.jar"),
        "REF_DIR": _resolve(ref_dir or os.environ.get("REF_DIR"), ROOT / "ref" / "b37.bref3"),
        "MAP_DIR": _resolve(map_dir or os.environ.get("MAP_DIR"), ROOT / "ref" / "maps"),
        "FASTA": _resolve(fasta or os.environ.get("FASTA"), ROOT / "ref" / "human_g1k_v37.fasta"),
    }
    missing_paths = [f"{name}={path}" for name, path in resolved.items() if not path.exists()]
    if missing_paths:
        raise RuntimeError("Chýbajú referenčné dáta pre imputáciu: " + "; ".join(missing_paths))
    if not any(resolved["REF_DIR"].glob("*.bref3")):
        raise RuntimeError(f"Referenčný priečinok neobsahuje bref3 súbory: {resolved['REF_DIR']}")

    work_dir.mkdir(parents=True, exist_ok=True)
    env = {
        **os.environ,
        **{name: str(path) for name, path in resolved.items()},
        "WORK_DIR": str(work_dir.resolve()),
        # The shell script historically called this DR2MIN, but extraction is
        # based on max(GP), which is the meaningful per-sample confidence.
        "DR2MIN": str(min_gp),
        "FLANK": str(int(flank)),
    }
    cmd = [
        bash, str(script), str(normalized), str(targets), str(status), str(out_prefix)
    ]
    proc = subprocess.run(
        cmd,
        cwd=str(ROOT),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log_path = Path(str(out_prefix) + "_imputation.log")
    log_path.write_text(proc.stdout or "", encoding="utf-8")
    if proc.returncode != 0:
        raise RuntimeError(
            f"Imputácia skončila s kódom {proc.returncode}. Detaily: {log_path}"
        )
    final = Path(str(out_prefix) + "_panel_final.csv")
    if not final.exists():
        raise RuntimeError(f"Imputácia nevytvorila očakávaný výstup: {final}")
    return final


def prepare_for_analyses(
    sample_path: str,
    output_dir: str,
    panel_path: Optional[str] = None,
    id_reference_path: Optional[str] = None,
    v5_template_path: Optional[str] = None,
    *,
    online_annotation: bool = False,
    run_imputation: bool = False,
    min_gp: float = 0.9,
    flank: int = 250_000,
    beagle_jar: Optional[str] = None,
    ref_dir: Optional[str] = None,
    map_dir: Optional[str] = None,
    fasta: Optional[str] = None,
) -> PreparationResult:
    """Normalize a sample once and produce sibling and predisposition outputs.

    ``v5_template_path`` is optional because the repository currently contains
    a v3/v4/v5 position union, not an exact 23andMe v5 manifest.  Without the
    template the sibling output is a v5-compatible GRCh37/plus-strand file that
    retains every measured marker.  With the template it is an exact marker-set
    scaffold and missing v5 markers are written as no-calls.
    """

    if not 0.0 <= float(min_gp) <= 1.0:
        raise ValueError("min_gp musí byť v intervale 0 až 1")
    if int(flank) < 0:
        raise ValueError("flank nesmie byť záporný")

    sample = Path(sample_path).expanduser().resolve()
    out_dir = Path(output_dir).expanduser().resolve()
    panel = _resolve(panel_path, DEFAULT_PANEL)
    id_reference = _resolve(id_reference_path, DEFAULT_ID_REFERENCE)
    for label, path in (("vzorka", sample), ("predispozičný panel", panel)):
        if not path.exists():
            raise FileNotFoundError(f"Chýba {label}: {path}")
    if not id_reference.exists():
        id_reference = None

    out_dir.mkdir(parents=True, exist_ok=True)
    stem = _safe_stem(str(sample))

    data = parse_file(str(sample))
    if data.build not in ("GRCh37", "b37", "hg19", "unknown", None):
        raise ValueError(
            f"Vzorka je označená ako {data.build}; pipeline vyžaduje GRCh37/hg19."
        )
    if data.source_format == "illumina_final_report":
        unsafe_strand = [w for w in data.warnings if "Treba strand-flip" in w]
        if unsafe_strand:
            raise ValueError(
                unsafe_strand[0]
                + " Pre spoločný sibling/predispozičný dataset treba export so stĺpcami "
                  "'Allele1 - Plus' a 'Allele2 - Plus' alebo validovanú manifestovú konverziu."
            )
        if id_reference:
            relabel_from_reference(data, str(id_reference))
        annotate_by_rsid(
            data,
            str(id_reference) if id_reference else None,
            online=online_annotation,
        )

    normalized = out_dir / f"{stem}.normalized_grch37_23andme.txt"
    write_23andme(data, str(normalized), drop_no_call=False, keep_non_rs=True)

    canonical_profile = "23andMe-v5-compatible"
    sibling = normalized
    sibling_markers = data.n_variants
    if v5_template_path:
        template = Path(v5_template_path).expanduser().resolve()
        if not template.exists():
            raise FileNotFoundError(f"Chýba šablóna 23andMe v5: {template}")
        scaffold, scaffold_stats = build_scaffold(data, str(template))
        sibling = out_dir / f"{stem}.sibling_23andme_v5.txt"
        write_23andme(scaffold, str(sibling), sort_alleles=False)
        sibling_markers = int(scaffold_stats["template_markers"])
        canonical_profile = "23andMe-v5-exact-template"

    panel_result = check_panel(
        data,
        str(panel),
        id_reference=str(id_reference) if id_reference else None,
        online=online_annotation,
    )
    status = out_dir / f"{stem}.predisposition_status.csv"
    targets = out_dir / f"{stem}.predisposition_missing.csv"
    panel_result.table.to_csv(status, index=False)
    panel_result.table[panel_result.table["status"] != "called"].to_csv(targets, index=False)

    out_prefix = out_dir / f"{stem}.predisposition"
    final_panel = out_dir / f"{stem}.predisposition_panel.csv"
    predisposition_dataset: Optional[Path] = None
    imputation_run = False
    rejected_low_conf = 0
    target_rows = pd.read_csv(targets, dtype=str).fillna("")
    imputable = target_rows[
        (target_rows.get("chromosome", "") != "") &
        (target_rows.get("position", "") != "")
    ] if len(target_rows) else target_rows

    if run_imputation and len(imputable):
        generated = _run_targeted_imputation(
            normalized=normalized,
            targets=targets,
            status=status,
            out_prefix=out_prefix,
            work_dir=out_dir / ".work" / stem,
            min_gp=float(min_gp),
            flank=int(flank),
            beagle_jar=beagle_jar,
            ref_dir=ref_dir,
            map_dir=map_dir,
            fasta=fasta,
        )
        rejected_low_conf = _filter_imputed_panel(generated, final_panel, float(min_gp))
        predisposition_dataset = out_dir / f"{stem}.with_predisposition_imputation.txt"
        merge_imputed_into_23andme(
            str(normalized), str(final_panel), str(predisposition_dataset), min_conf=float(min_gp)
        )
        imputation_run = True
    else:
        merge_panel(str(status)).to_csv(final_panel, index=False)

    counts = _panel_counts(final_panel)
    report_path = out_dir / f"{stem}.preparation_report.json"
    result = PreparationResult(
        sample=str(sample),
        source_format=data.source_format,
        canonical_profile=canonical_profile,
        normalized_dataset=str(normalized),
        sibling_dataset=str(sibling),
        predisposition_panel=str(final_panel),
        predisposition_dataset=str(predisposition_dataset) if predisposition_dataset else None,
        panel_status=str(status),
        panel_missing=str(targets),
        report=str(report_path),
        measured_variants=int(data.n_variants),
        measured_called=int(data.n_called),
        sibling_markers=int(sibling_markers),
        panel_markers=int(panel_result.n_panel),
        panel_measured=int(counts["measured"]),
        panel_imputed=int(counts["imputed"]),
        panel_rejected_low_conf=int(rejected_low_conf),
        panel_missing_after=int(counts["missing"]),
        imputation_run=imputation_run,
    )
    report = asdict(result)
    report["reference_build"] = "GRCh37"
    report["sibling_contains_imputed_genotypes"] = False
    report["warnings"] = list(data.warnings)
    if not v5_template_path:
        report["warnings"].append(
            "Nebola dodaná presná šablóna 23andMe v5; sibling výstup používa "
            "v5-kompatibilný formát a zachováva všetky namerané markery."
        )
    if not run_imputation and len(imputable):
        report["warnings"].append(
            f"Imputácia nebola spustená; {len(imputable)} panelových lokusov má známu pozíciu."
        )
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return result

