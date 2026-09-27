"""
Unit and Integration Tests for cas12a_xai.
"""

import unittest
import numpy as np
import pandas as pd

from cas12a_xai.features import (
    calc_gc,
    calc_tm,
    calc_stacking_dg,
    decompose_target_sequence,
    extract_features_single,
    extract_features_batch,
    FEATURE_NAMES
)
from cas12a_xai.predictor import Cas12aPredictor, predict_efficiency, explain_guide


class TestCas12aFeatureExtraction(unittest.TestCase):

    def setUp(self):
        # Canonical test sequences
        self.seq_34bp = "ACGGTTTCGACCGTTAGGCTACGATCGGATCGCC"  # 4bp flank + TTTC + 23bp + 3bp flank = 34bp
        self.seq_27bp = "TTTCGACCGTTAGGCTACGATCGGATC"          # TTTC + 23bp = 27bp
        self.seq_23bp = "GACCGTTAGGCTACGATCGGATC"              # 23bp spacer alone
        self.seq_poly_t = "TTTCAAAATTTTCCCGGGAAACCCGAA"        # 4bp PAM + 23bp spacer containing TTTT = 27bp

    def test_calc_gc(self):
        self.assertAlmostEqual(calc_gc("GGCC"), 1.0)
        self.assertAlmostEqual(calc_gc("AATT"), 0.0)
        self.assertAlmostEqual(calc_gc("GCAT"), 0.5)

    def test_calc_tm(self):
        # Wallace rule (<14 nt): (wA+xT)*2 + (yG+zC)*4
        # "ATGC" has 2 AT and 2 GC -> 2*2 + 2*4 = 12
        self.assertEqual(calc_tm("ATGC"), 12.0)
        # Marmur rule (>=14 nt): 64.9 + 41*(GC - 16.4)/len
        tm_long = calc_tm("ATGCATGCATGCATGC")
        self.assertTrue(30.0 < tm_long < 70.0)

    def test_calc_stacking_dg(self):
        # "AA" is -1.00 kcal/mol
        self.assertAlmostEqual(calc_stacking_dg("AA"), -1.00)
        # "AAAA" has 3 AA steps -> -3.00 kcal/mol
        self.assertAlmostEqual(calc_stacking_dg("AAAA"), -3.00)
        # "CG" is -2.17 kcal/mol
        self.assertAlmostEqual(calc_stacking_dg("CG"), -2.17)

    def test_sequence_decomposition(self):
        d34 = decompose_target_sequence(self.seq_34bp)
        self.assertEqual(len(d34['upstream_4bp']), 4)
        self.assertEqual(d34['pam_4bp'], "TTTC")
        self.assertEqual(len(d34['spacer_23bp']), 23)
        self.assertEqual(len(d34['seed_8bp']), 8)
        self.assertEqual(len(d34['trunk_8bp']), 8)
        self.assertEqual(len(d34['distal_7bp']), 7)

        d27 = decompose_target_sequence(self.seq_27bp)
        self.assertEqual(d27['pam_4bp'], "TTTC")
        self.assertEqual(len(d27['spacer_23bp']), 23)

    def test_feature_vector_dimension(self):
        feats = extract_features_single(self.seq_34bp)
        self.assertEqual(len(feats), len(FEATURE_NAMES))
        for feat_name in FEATURE_NAMES:
            self.assertIn(feat_name, feats)

    def test_poly_t_terminator_detection(self):
        feats_normal = extract_features_single(self.seq_27bp)
        self.assertEqual(feats_normal['poly_t_terminator'], 0)

        feats_poly_t = extract_features_single(self.seq_poly_t)
        self.assertEqual(feats_poly_t['poly_t_terminator'], 1)

    def test_batch_extraction(self):
        df = extract_features_batch([self.seq_34bp, self.seq_27bp])
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 2)
        self.assertEqual(list(df.columns), FEATURE_NAMES)


class TestCas12aPredictor(unittest.TestCase):

    def setUp(self):
        self.predictor = Cas12aPredictor()
        self.active_seq = "TTTCGACCGTTAGGCTACGATCGGATC"
        self.inactive_seq = "TTTCAAAATTTTCCCGGGAAACCCGAA"

    def test_prediction_output_range(self):
        pred = self.predictor.predict(self.active_seq)
        self.assertIsInstance(pred, np.ndarray)
        self.assertTrue(0.0 <= pred[0] <= 100.0)

    def test_poly_t_repression(self):
        pred_active = self.predictor.predict(self.active_seq)[0]
        pred_inactive = self.predictor.predict(self.inactive_seq)[0]
        self.assertGreater(pred_active, pred_inactive)

    def test_explain_guide(self):
        explanation = self.predictor.explain(self.inactive_seq)
        self.assertIn("target_sequence", explanation)
        self.assertIn("predicted_indel_percent", explanation)
        self.assertIn("biophysical_metrics", explanation)
        self.assertIn("flags", explanation)
        self.assertIn("design_recommendations", explanation)
        # Should flag poly-T terminator and override activity tier
        flag_str = " ".join(explanation['flags'])
        self.assertIn("Poly-T", flag_str)
        self.assertIn("Construct-Risk Critical", explanation['activity_tier'])

    def test_boundary_poly_t_escape(self):
        """
        Verify biological boundary escape:
        If TTTT spans positions 19-22 (across the 20-nt guide / 3-bp target context junction),
        the 20-nt guide transcript only contains 2 Ts ('TT') and escapes premature termination.
        """
        boundary_seq = "TTTCGACCGTTAGGCTACGATCTTTTC"
        explanation = self.predictor.explain(boundary_seq)
        self.assertEqual(explanation["parsed_domains"]["has_guide20_polyt"], 0)
        self.assertEqual(explanation["parsed_domains"]["has_boundary_polyt"], 1)
        self.assertIn("High Cleavage", explanation["activity_tier"])
        flag_str = " ".join(explanation["flags"])
        self.assertIn("Boundary Poly-T", flag_str)

    def test_cpu_inference_latency_benchmark(self):
        """
        Verify that inference latency meets high-throughput production standards:
        - End-to-end predict() for 1,000 guides executes in < 3,000 ms (< 3 ms/guide) on CPU.
        - If a trained booster model is loaded, also test booster predict latency (< 100 ms).
        """
        import time
        guides_1k = [self.active_seq] * 1000

        # Full end-to-end predict() latency (feature extraction + scoring)
        start_e2e = time.perf_counter()
        preds_e2e = self.predictor.predict(guides_1k)
        e2e_elapsed_ms = (time.perf_counter() - start_e2e) * 1000.0
        self.assertEqual(len(preds_e2e), 1000)
        self.assertLess(e2e_elapsed_ms, 3000.0, f"E2E predict latency too slow: {e2e_elapsed_ms:.2f} ms")

        # If a trained booster model is loaded, verify pure booster latency
        if self.predictor.model is not None:
            X_df = extract_features_batch(guides_1k)
            start_inf = time.perf_counter()
            preds_booster = self.predictor.model.predict(X_df.values)
            inf_elapsed_ms = (time.perf_counter() - start_inf) * 1000.0
            self.assertEqual(len(preds_booster), 1000)
            self.assertLess(inf_elapsed_ms, 100.0, f"Pure inference latency too slow: {inf_elapsed_ms:.2f} ms")


if __name__ == "__main__":
    unittest.main()
