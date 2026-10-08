"""Security and resource guards shared by the local and hosted application."""

from __future__ import annotations

import base64
import os
import secrets
import threading
import time
from pathlib import Path
from typing import Optional

from fastapi import HTTPException, UploadFile


def env_flag(name: str, default: bool = False) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


PUBLIC_MODE = env_flag("PUBLIC_MODE")
MAX_UPLOAD_MB = int(os.environ.get("MAX_UPLOAD_MB", "100" if PUBLIC_MODE else "2048"))
MAX_REQUEST_MB = int(os.environ.get("MAX_REQUEST_MB", "250" if PUBLIC_MODE else "4096"))
RETENTION_HOURS = float(os.environ.get("RUNTIME_RETENTION_HOURS", "6" if PUBLIC_MODE else "168"))
MAX_BATCH_SAMPLES = int(os.environ.get("MAX_BATCH_SAMPLES", "5" if PUBLIC_MODE else "100"))

# Beagle and the batch runner share scratch/cache locations.  One process-wide
# lock avoids cross-job contamination and also caps public compute consumption.
JOB_LOCK = threading.Lock()


def validate_public_config() -> None:
    """Fail closed when a public deployment has no authentication secret."""

    if not PUBLIC_MODE:
        return
    missing = [name for name in ("APP_USERNAME", "APP_PASSWORD") if not os.environ.get(name)]
    if missing:
        raise RuntimeError(
            "PUBLIC_MODE=1 vyžaduje Railway secrets: " + ", ".join(missing)
        )
    demo_values = [os.environ.get("DEMO_USERNAME"), os.environ.get("DEMO_PASSWORD")]
    if any(demo_values) and not all(demo_values):
        raise RuntimeError(
            "DEMO účet vyžaduje obe premenné: DEMO_USERNAME a DEMO_PASSWORD"
        )


def basic_auth_role(authorization: Optional[str]) -> Optional[str]:
    """Return ``owner`` or ``demo`` after constant-time Basic Auth validation."""

    if not PUBLIC_MODE:
        return "owner"
    if not authorization or not authorization.startswith("Basic "):
        return None
    try:
        raw = base64.b64decode(authorization[6:], validate=True).decode("utf-8")
        username, password = raw.split(":", 1)
    except (ValueError, UnicodeDecodeError):
        return None
    expected_user = os.environ.get("APP_USERNAME", "")
    expected_password = os.environ.get("APP_PASSWORD", "")
    if secrets.compare_digest(username, expected_user) and secrets.compare_digest(
        password, expected_password
    ):
        return "owner"
    demo_user = os.environ.get("DEMO_USERNAME", "")
    demo_password = os.environ.get("DEMO_PASSWORD", "")
    if demo_user and demo_password:
        if secrets.compare_digest(username, demo_user) and secrets.compare_digest(
            password, demo_password
        ):
            return "demo"
    return None


def basic_credentials_valid(authorization: Optional[str]) -> bool:
    """Backward-compatible boolean wrapper used by callers and tests."""

    return basic_auth_role(authorization) is not None


def save_upload_limited(
    upload: UploadFile,
    destination: Path,
    *,
    max_bytes: Optional[int] = None,
) -> Path:
    """Stream an upload to disk and remove partial data when it exceeds the limit."""

    limit = int(max_bytes if max_bytes is not None else MAX_UPLOAD_MB * 1024 * 1024)
    destination.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    try:
        with destination.open("wb") as handle:
            while True:
                chunk = upload.file.read(1024 * 1024)
                if not chunk:
                    break
                written += len(chunk)
                if written > limit:
                    raise HTTPException(
                        status_code=413,
                        detail=f"Súbor prekročil povolený limit {MAX_UPLOAD_MB} MB.",
                    )
                handle.write(chunk)
    except Exception:
        destination.unlink(missing_ok=True)
        raise
    return destination


def prune_runtime_files(directory: Path, *, now: Optional[float] = None) -> int:
    """Delete expired runtime files under one explicitly supplied work directory."""

    if not directory.exists():
        return 0
    cutoff = (time.time() if now is None else now) - RETENTION_HOURS * 3600
    removed = 0
    files = [path for path in directory.rglob("*") if path.is_file() and not path.is_symlink()]
    for path in files:
        try:
            if path.stat().st_mtime < cutoff:
                path.unlink()
                removed += 1
        except FileNotFoundError:
            pass
    directories = sorted(
        (path for path in directory.rglob("*") if path.is_dir() and not path.is_symlink()),
        key=lambda path: len(path.parts),
        reverse=True,
    )
    for path in directories:
        try:
            path.rmdir()
        except OSError:
            pass
    return removed
