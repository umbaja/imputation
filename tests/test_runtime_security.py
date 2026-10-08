import base64
import io
import os
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock

from fastapi import HTTPException, UploadFile

from app import runtime_security


class RuntimeSecurityTests(unittest.TestCase):
    def test_save_upload_under_limit(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "sample.txt"
            upload = UploadFile(filename="sample.txt", file=io.BytesIO(b"ACGT"))
            runtime_security.save_upload_limited(upload, destination, max_bytes=4)
            self.assertEqual(destination.read_bytes(), b"ACGT")

    def test_save_upload_over_limit_removes_partial_file(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "sample.txt"
            upload = UploadFile(filename="sample.txt", file=io.BytesIO(b"ACGTA"))
            with self.assertRaises(HTTPException) as error:
                runtime_security.save_upload_limited(upload, destination, max_bytes=4)
            self.assertEqual(error.exception.status_code, 413)
            self.assertFalse(destination.exists())

    def test_prune_removes_only_expired_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            old = root / "old.txt"
            fresh = root / "fresh.txt"
            old.write_text("old", encoding="utf-8")
            fresh.write_text("fresh", encoding="utf-8")
            now = time.time()
            os.utime(old, (now - 7200, now - 7200))
            with mock.patch.object(runtime_security, "RETENTION_HOURS", 1):
                removed = runtime_security.prune_runtime_files(root, now=now)
            self.assertEqual(removed, 1)
            self.assertFalse(old.exists())
            self.assertTrue(fresh.exists())

    def test_basic_auth_in_public_mode(self):
        token = base64.b64encode(b"tester:long-secret").decode("ascii")
        environment = {"APP_USERNAME": "tester", "APP_PASSWORD": "long-secret"}
        with mock.patch.object(runtime_security, "PUBLIC_MODE", True), \
                mock.patch.dict(os.environ, environment, clear=False):
            self.assertTrue(runtime_security.basic_credentials_valid(f"Basic {token}"))
            self.assertFalse(runtime_security.basic_credentials_valid("Basic invalid"))
            self.assertFalse(runtime_security.basic_credentials_valid(None))

    def test_public_mode_fails_closed_without_credentials(self):
        with mock.patch.object(runtime_security, "PUBLIC_MODE", True), \
                mock.patch.dict(os.environ, {}, clear=True):
            with self.assertRaises(RuntimeError):
                runtime_security.validate_public_config()


if __name__ == "__main__":
    unittest.main()
