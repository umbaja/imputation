import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class RailwayConfigTests(unittest.TestCase):
    def test_start_script_binds_to_railway_port_and_can_bootstrap_reference(self):
        script = (ROOT / "deploy" / "start.sh").read_text(encoding="utf-8")
        self.assertIn('${PORT:-8000}', script)
        self.assertIn("BOOTSTRAP_REFERENCE", script)
        self.assertIn("setup_reference.sh", script)

    def test_dockerfile_runs_start_script(self):
        dockerfile = (ROOT / "Dockerfile").read_text(encoding="utf-8")
        self.assertIn('CMD ["bash", "deploy/start.sh"]', dockerfile)

    def test_application_exposes_healthcheck(self):
        main = (ROOT / "app" / "main.py").read_text(encoding="utf-8")
        self.assertIn('@app.get("/healthz")', main)


if __name__ == "__main__":
    unittest.main()
