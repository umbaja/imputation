import csv
import re
import unittest
from pathlib import Path, PurePosixPath


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "reference_manifest.tsv"


class ReferenceManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with MANIFEST.open(encoding="utf-8", newline="") as handle:
            cls.rows = list(csv.DictReader(handle, delimiter="\t"))

    def test_schema_and_unique_paths(self):
        self.assertEqual(
            set(self.rows[0]),
            {"kind", "chromosome", "path", "size_bytes", "sha256", "url"},
        )
        paths = [row["path"] for row in self.rows]
        self.assertEqual(len(paths), len(set(paths)))

    def test_expected_components_are_pinned(self):
        kinds = [row["kind"] for row in self.rows]
        self.assertEqual(kinds.count("beagle"), 1)
        self.assertEqual(kinds.count("maps_archive"), 1)
        self.assertEqual(kinds.count("fasta_archive"), 1)
        self.assertEqual(kinds.count("fasta"), 1)

        bref_chromosomes = {
            row["chromosome"] for row in self.rows if row["kind"] == "bref3"
        }
        self.assertEqual(bref_chromosomes, {*(str(i) for i in range(1, 23)), "X"})

    def test_checksums_sizes_urls_and_paths(self):
        for row in self.rows:
            with self.subTest(path=row["path"]):
                self.assertRegex(row["sha256"], re.compile(r"^[0-9a-f]{64}$"))
                self.assertGreater(int(row["size_bytes"]), 0)

                path = PurePosixPath(row["path"])
                self.assertFalse(path.is_absolute())
                self.assertNotIn("..", path.parts)

                if row["kind"] == "fasta":
                    self.assertEqual(row["url"], "-")
                else:
                    self.assertTrue(row["url"].startswith("https://"))

    def test_reference_directory_stays_ignored(self):
        ignore_lines = {
            line.strip() for line in (ROOT / ".gitignore").read_text().splitlines()
        }
        self.assertIn("ref/", ignore_lines)


if __name__ == "__main__":
    unittest.main()
