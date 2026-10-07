from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from dna_tools.analysis_prep import _filter_imputed_panel, prepare_for_analyses


class AnalysisPreparationTests(unittest.TestCase):
    def test_low_confidence_imputation_is_not_reported_as_a_call(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.csv"
            destination = root / "filtered.csv"
            pd.DataFrame([
                {"rsid": "rs1", "chromosome": "1", "position": "10", "genotype": "AA",
                 "source": "measured", "confidence": "", "dr2": ""},
                {"rsid": "rs2", "chromosome": "1", "position": "20", "genotype": "AG",
                 "source": "imputed", "confidence": "0.95", "dr2": "0.8"},
                {"rsid": "rs3", "chromosome": "1", "position": "30", "genotype": "GG",
                 "source": "imputed", "confidence": "0.60", "dr2": "0.5"},
                {"rsid": "rs4", "chromosome": "1", "position": "40", "genotype": "CT",
                 "source": "imputed", "confidence": "", "dr2": ""},
            ]).to_csv(source, index=False)

            rejected = _filter_imputed_panel(source, destination, min_gp=0.9)
            result = pd.read_csv(destination, dtype=str, keep_default_na=False)

            self.assertEqual(rejected, 2)
            self.assertEqual(result.loc[result.rsid == "rs2", "source"].iloc[0], "imputed")
            self.assertEqual(result.loc[result.rsid == "rs3", "source"].iloc[0], "missing")
            self.assertEqual(result.loc[result.rsid == "rs3", "genotype"].iloc[0], "--")
            self.assertEqual(result.loc[result.rsid == "rs4", "source"].iloc[0], "missing")

    def test_one_normalization_creates_separate_analysis_products(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample = root / "sample.txt"
            sample.write_text(
                "# rsid\tchromosome\tposition\tgenotype\n"
                "rs1\t1\t10\tAA\n"
                "rs2\t1\t20\tCT\n",
                encoding="utf-8",
            )
            panel = root / "panel.csv"
            panel.write_text(
                "rsid,chromosome,position\n"
                "rs1,1,10\n"
                "rs2,1,20\n"
                "rs3,1,30\n",
                encoding="utf-8",
            )
            output = root / "output"

            result = prepare_for_analyses(
                str(sample),
                str(output),
                panel_path=str(panel),
                id_reference_path=str(root / "missing-reference.tsv.gz"),
            )

            self.assertEqual(result.panel_markers, 3)
            self.assertEqual(result.panel_measured, 2)
            self.assertEqual(result.panel_imputed, 0)
            self.assertEqual(result.panel_missing_after, 1)
            self.assertEqual(result.panel_rejected_low_conf, 0)
            self.assertEqual(result.sibling_dataset, result.normalized_dataset)
            report = json.loads(Path(result.report).read_text(encoding="utf-8"))
            self.assertFalse(report["sibling_contains_imputed_genotypes"])

    def test_top_strand_illumina_export_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sample = root / "illumina.txt"
            sample.write_text(
                "[Header]\nContent\tTest\n[Data]\n"
                "SNP Name\tChr\tPosition\tAllele1 - Top\tAllele2 - Top\n"
                "rs1\t1\t10\tA\tG\n",
                encoding="utf-8",
            )
            panel = root / "panel.csv"
            panel.write_text("rsid,chromosome,position\nrs1,1,10\n", encoding="utf-8")

            with self.assertRaisesRegex(ValueError, "strand-flip"):
                prepare_for_analyses(
                    str(sample),
                    str(root / "output"),
                    panel_path=str(panel),
                    id_reference_path=str(root / "missing-reference.tsv.gz"),
                )


if __name__ == "__main__":
    unittest.main()
