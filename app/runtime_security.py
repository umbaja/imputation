"""Security and resource guards shared by the local and hosted application."""

from __future__ import annotations

import base64
import json
import os
import secrets
import threading
import time
from datetime import datetime
from pathlib import Path
from typing import Optional
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

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
TEST_DAILY_LIMIT = int(os.environ.get("TEST_DAILY_LIMIT", "5"))
TEST_QUOTA_TIMEZONE = os.environ.get("TEST_QUOTA_TIMEZONE", "Europe/Bratislava")
TEST_QUOTA_FILE = Path(
    os.environ.get(
        "TEST_QUOTA_FILE",
        "/app/ref/test_account_quota.json" if PUBLIC_MODE else "work/test_account_quota.json",
    )
)

# Beagle and the batch runner share scratch/cache locations.  One process-wide
# lock avoids cross-job contamination and also caps public compute consumption.
JOB_LOCK = threading.Lock()
TEST_QUOTA_LOCK = threading.RLock()


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
    test_values = [os.environ.get("TEST_USERNAME"), os.environ.get("TEST_PASSWORD")]
    if any(test_values) and not all(test_values):
        raise RuntimeError(
            "Testovací účet vyžaduje obe premenné: TEST_USERNAME a TEST_PASSWORD"
        )
    if TEST_DAILY_LIMIT < 1:
        raise RuntimeError("TEST_DAILY_LIMIT musí byť kladné celé číslo")


def basic_auth_role(authorization: Optional[str]) -> Optional[str]:
    """Return an account role after constant-time Basic Auth validation."""

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
    test_user = os.environ.get("TEST_USERNAME", "")
    test_password = os.environ.get("TEST_PASSWORD", "")
    if test_user and test_password:
        if secrets.compare_digest(username, test_user) and secrets.compare_digest(
            password, test_password
        ):
            return "test"
    return None


def _quota_day() -> str:
    try:
        timezone = ZoneInfo(TEST_QUOTA_TIMEZONE)
    except ZoneInfoNotFoundError as exc:
        raise RuntimeError(
            f"Neznáme časové pásmo TEST_QUOTA_TIMEZONE={TEST_QUOTA_TIMEZONE!r}"
        ) from exc
    return datetime.now(timezone).date().isoformat()


def _read_test_quota_unlocked(day: str) -> dict:
    try:
        payload = json.loads(TEST_QUOTA_FILE.read_text(encoding="utf-8"))
    except (FileNotFoundError, OSError, ValueError, TypeError):
        payload = {}
    if payload.get("day") != day:
        return {"day": day, "used": 0}
    try:
        used = max(0, int(payload.get("used", 0)))
    except (TypeError, ValueError):
        used = 0
    return {"day": day, "used": used}


def _write_test_quota_unlocked(payload: dict) -> None:
    TEST_QUOTA_FILE.parent.mkdir(parents=True, exist_ok=True)
    temporary = TEST_QUOTA_FILE.with_suffix(TEST_QUOTA_FILE.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(TEST_QUOTA_FILE)


def test_quota_status() -> dict:
    """Return the shared test account's persistent quota for the current day."""

    day = _quota_day()
    with TEST_QUOTA_LOCK:
        payload = _read_test_quota_unlocked(day)
    used = min(payload["used"], TEST_DAILY_LIMIT)
    return {
        "day": day,
        "limit": TEST_DAILY_LIMIT,
        "used": used,
        "remaining": max(0, TEST_DAILY_LIMIT - used),
        "timezone": TEST_QUOTA_TIMEZONE,
    }


def reserve_test_quota(role: str, count: int) -> dict:
    """Atomically reserve genotype slots for the shared test role."""

    if role != "test":
        return {}
    if count < 1:
        raise HTTPException(status_code=400, detail="Počet genotypov musí byť kladný.")
    day = _quota_day()
    with TEST_QUOTA_LOCK:
        payload = _read_test_quota_unlocked(day)
        remaining = max(0, TEST_DAILY_LIMIT - payload["used"])
        if count > remaining:
            raise HTTPException(
                status_code=429,
                detail=(
                    f"Testovací účet má spoločný limit {TEST_DAILY_LIMIT} genotypov denne; "
                    f"dnes zostáva {remaining}. Limit sa obnoví o polnoci "
                    f"({TEST_QUOTA_TIMEZONE})."
                ),
            )
        payload["used"] += count
        _write_test_quota_unlocked(payload)
    return test_quota_status()


def refund_test_quota(role: str, count: int) -> None:
    """Return a reservation when a batch could not be accepted for processing."""

    if role != "test" or count < 1:
        return
    day = _quota_day()
    with TEST_QUOTA_LOCK:
        payload = _read_test_quota_unlocked(day)
        payload["used"] = max(0, payload["used"] - count)
        _write_test_quota_unlocked(payload)


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
