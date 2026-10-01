"""
cas12a_xai: Explainable Machine Learning for CRISPR-Cas12a On-Target Cleavage Efficiency.

Developed by Minh Tran (Independent Researcher, Perth, Western Australia, Australia).
Integrates SantaLucia nearest-neighbor thermodynamics, seed-distal polarity gradients,
domain-segmented feature engineering, and TreeSHAP feature attributions.
"""

from .features import (
    extract_features_single,
    extract_features_batch,
    calc_gc,
    calc_tm,
    calc_stacking_dg,
    NN_DELTA_G,
    FEATURE_NAMES
)
from .predictor import Cas12aPredictor, predict_efficiency, explain_guide

__version__ = "0.1.0"
__author__ = "Minh Tran"
__license__ = "MIT"

__all__ = [
    "extract_features_single",
    "extract_features_batch",
    "calc_gc",
    "calc_tm",
    "calc_stacking_dg",
    "NN_DELTA_G",
    "FEATURE_NAMES",
    "Cas12aPredictor",
    "predict_efficiency",
    "explain_guide",
]
