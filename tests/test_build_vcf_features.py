from __future__ import annotations

import gzip
import tempfile
import unittest
from pathlib import Path

from src.features.build_vcf_features import build_features


HEADER = (
    "chrom\tpos\tid\tref\talt\tqual\tfilter\tac\taf\tan\tinfo_dp\tfs\tmq\tqd\t"
    "mq_rank_sum\tread_pos_rank_sum\tbase_q_rank_sum\tclipping_rank_sum\tsor\t"
    "excess_het\tinbreeding_coeff\tgt\tad\tsample_dp\tgq\tpl\tpgt\tpid\n"
)


class BuildVcfFeaturesTest(unittest.TestCase):
    def test_derives_depth_vaf_and_variant_flags(self) -> None:
        row = (
            "1\t100\trs123\tA\tT\t120\tPASS\t1\t0.5\t2\t30\t1\t60\t10\t0\t0\t0\t0\t1\t0\t0\t"
            "0/1\t24,6\t30\t80\t.\t0|1\t100_A_T\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "features.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(HEADER)
                handle.write(row)

            result = build_features(path).collect().row(0, named=True)

        self.assertEqual(result["ref_depth"], 24)
        self.assertEqual(result["alt_depth"], 6)
        self.assertAlmostEqual(result["vaf"], 0.2)
        self.assertTrue(result["is_pass"])
        self.assertTrue(result["is_snp"])
        self.assertTrue(result["is_het_ref_alt"])
        self.assertTrue(result["is_dbsnp"])
        self.assertTrue(result["has_local_phase"])
        self.assertEqual(result["submission_chrom"], "chr1")
        self.assertEqual(result["variant_key"], "chr1:100:A:T")


if __name__ == "__main__":
    unittest.main()
