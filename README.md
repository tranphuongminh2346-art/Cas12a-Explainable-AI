# cas12a-xai: Biophysically-Guided Explainable Machine Learning for CRISPR-Cas12a Cleavage Efficiency

[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Status: In Review](https://img.shields.io/badge/Manuscript-Targeting%20Q1%2FQ2%20Bioinformatics-green.svg)](#citation)
[![Architecture: LightGBM + TreeSHAP](https://img.shields.io/badge/Architecture-LightGBM%20%2B%20TreeSHAP-orange.svg)](#overview)
[![Tests](https://img.shields.io/badge/Tests-Passing-brightgreen.svg)](#testing)

> **Author**: Minh Tran  
> *School of Biomedical Sciences, University of Western Australia / Independent Researcher*  
> *Correspondence: minh.tran@research.uwa.edu.au*  

---

## Overview

**cas12a-xai** is an open-source, biophysically-grounded machine learning framework and design toolkit for predicting and optimizing CRISPR-Cas12a (Cpf1) guide RNA on-target cleavage efficiency. 

Unlike opaque deep neural networks that act as "black boxes", **cas12a-xai** models crRNA hybridization using domain-segmented SantaLucia (1998) nearest-neighbor thermodynamics, kinetic R-loop polarity gradients, and position-specific nucleotide encodings. By pairing gradient boosted decision trees (LightGBM) with exact Shapley additive explanations (TreeSHAP), this toolkit delivers both **high predictive fidelity** and **transparent, mechanistic biological insights** into the molecular determinants of Cas12a cutting.

```
                    CRISPR-Cas12a TARGET SITE (34 bp)
 ┌─────────┬──────────────┬─────────────────────────┬──────────────┐
 │ 5' Flank│  PAM (4 bp)  │     SPACER (23 bp)      │ 3' Flank     │
 │ (-8..-5)│   5'-TTTV    │ Seed(1-8) Trunk(9-16) Distal(17-23)    │ (24..26)     │
 └─────────┴──────────────┴─────────────────────────┴──────────────┘
        │          │                   │
        ▼          ▼                   ▼
 [Upstream GC] [PAM PI Cleft] [Nearest-Neighbor ΔG°] ───► [LightGBM Ensemble]
                              [Seed-Distal ΔTm Grad]             │
                              [Poly-T Terminator]                ▼
                                                       [Predicted Indel %]
                                                       [TreeSHAP Attribution]
```

---

## Key Biophysical Discoveries

1. **Catastrophic Poly-T Repression (-35.2% Indel)**: Internal `TTTT` tracts act as intrinsic RNA Polymerase III transcription terminators, truncating nascent crRNA transcripts prior to ribonucleoprotein assembly.
2. **Duplex Base-Stacking Stabilization (+15.3% Indel)**: Favorable nearest-neighbor stacking free energies ($\Delta G^\circ_{\text{spacer}} < -32\text{ kcal/mol}$) promote stable R-loop formation and non-target strand displacement.
3. **Seed-to-Distal Kinetic Polarity Gradient (+8.4% Indel)**: A positive melting temperature gradient ($\Delta T_m = T_{m, \text{seed}} - T_{m, \text{distal}} > 0$) ensures low activation energy for nucleation while avoiding hyper-stable distal duplexes that stall the RuvC catalytic release.
4. **Epistatic PAM $\times$ Seed Architecture**: TTTC PAM allows flexible seed initiation, whereas suboptimal TTTG PAM demands rigid base stacking ($\Delta G^\circ_{\text{seed}} < -9.5\text{ kcal/mol}$) to achieve productive cleavage.
5. **Position 1 Thymine Repression (-10.1% Indel)**: A Thymine directly adjacent to the PAM severely destabilizes the initial unwinding bubble.

---

## Benchmark Highlights (Out-of-Fold Cross-Validation)

Evaluated across **11,365 canonical targets** (5'-TTTV-3' PAM) from human high-throughput screening data (Kim et al., 2018; 10-Fold Stratified Cross-Validation):

| Cohort | Sample Size ($n$) | Spearman $\rho$ | Pearson $r$ | ROC-AUC ($\ge 20\%$) | MAE (%) | RMSE (%) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **TTTA PAM** | 1,400 | **0.530** | **0.576** | 0.812 | 20.53 | 25.42 |
| **TTTC PAM** | 5,068 | **0.518** | **0.562** | 0.807 | 19.85 | 24.58 |
| **TTTG PAM** | 4,897 | **0.521** | **0.544** | 0.790 | 20.33 | 25.02 |
| **Combined Canonical** | **11,365** | **0.525** | **0.558** | **0.801** | **20.14** | **24.87** |

*Inference latency is hardware-dependent, typically hundreds of milliseconds per 1,000 guides (~150–600 ms, < 0.6 ms/guide) on standard consumer CPUs (no GPU required).*

---

## Installation

### Prerequisites
- Python 3.8, 3.9, 3.10, or 3.11
- Git

### Install via pip (Development / Editable Mode)
```bash
git clone https://github.com/minhtran-bio/cas12a-xai.git
cd cas12a-xai
pip install -e .
```

### Install Dependencies Directly
```bash
pip install -r requirements.txt
```

---

## Quickstart & Usage

### 1. Command-Line Interface (CLI)

#### Explain and Critique a Single Guide RNA
```bash
cas12a-xai --seq TTTCGACCGTTAGGCTACGATCGGATC --explain
```
**Example Output Report:**
```text
======================================================================
CRISPR-Cas12a Explainable AI Design Report
======================================================================
Target Sequence:       TTTCGACCGTTAGGCTACGATCGGATC
Predicted Indel Eff:   92.45% (High Cleavage (>=20%))
PAM Motif:             TTTC
Spacer (23 bp):        GACCGTTAGGCTACGATCGGATC
----------------------------------------------------------------------
Biophysical Properties:
  - spacer_gc_percent                         : 56.5
  - seed_gc_percent                           : 50.0
  - distal_gc_percent                         : 57.1
  - seed_tm_celsius                           : 24.0
  - distal_tm_celsius                         : 22.0
  - tm_polarity_gradient                      : 2.0
  - spacer_stacking_free_energy_kcal_per_mol  : -31.57
----------------------------------------------------------------------
Actionable Design Recommendations:
  [+] Optimal PAM subclass: TTTC maximizes catalytic cleft accommodation (+4.2% predicted gain).
  [+] Favorable seed initiation: Base 'G' at position 1 stabilizes the opening bubble.
  [+] Favorable kinetic polarity gradient (Delta_Tm = +2.0°C): Seed is thermodynamically more stable than distal region, facilitating rapid forward R-loop zippering.
======================================================================
```

#### Detect Premature Poly-T Termination
```bash
cas12a-xai --seq TTTCAAAATTTTCCCGGGAAACCCGAA --explain
```
```text
Target Sequence:       TTTCAAAATTTTCCCGGGAAACCCGAA
Predicted Indel Eff:   38.46% (Construct-Risk Critical (Premature Pol III Termination))
PAM Motif:             TTTC
Spacer (23 bp):        AAAATTTTCCCGGGAAACCCGAA
----------------------------------------------------------------------
Attention Flags:
  [!] CRITICAL: Premature Poly-T Pol III terminator detected in 20-nt guide transcript ('TTTT')
----------------------------------------------------------------------
Actionable Design Recommendations:
  [+] Severe expression ablation: nascent crRNA will be prematurely truncated by U6/7SK Pol III. Introduce synonymous mutations if targeting coding sequences, or shift guide offset.
======================================================================
```

#### High-Throughput Batch Screening (CSV)
```bash
cas12a-xai --file library_candidates.csv --column target_seq --output scored_library.csv
```

---

### 2. Python API

```python
from cas12a_xai import predict_efficiency, explain_guide, extract_features_single

# 1. Rapid Indel Cleavage Prediction (%)
sequence = "TTTCGACCGTTAGGCTACGATCGGATC"
indel_eff = predict_efficiency(sequence)
print(f"Predicted Indel Frequency: {indel_eff:.2f}%")

# 2. Comprehensive Explainability & Biophysical Breakdown
report = explain_guide(sequence)
print(f"Activity Tier: {report['activity_tier']}")
print(f"Polarity Gradient Delta_Tm: {report['biophysical_metrics']['tm_polarity_gradient']} °C")
for rec in report['design_recommendations']:
    print(f" -> {rec}")

# 3. Extract Raw Biophysical Feature Dictionary
features = extract_features_single(sequence)
print(f"Base Stacking Free Energy: {features['stacking_dg_spacer']} kcal/mol")

# 4. Custom GBDT Booster Inference (Optional)
# from cas12a_xai import Cas12aPredictor
# predictor = Cas12aPredictor(model_path="path/to/trained_booster.txt")
# custom_preds = predictor.predict([sequence])
```

> **Dual-Engine Architecture & Explainability**:
> 1. **Calibrated Biophysical Scoring Engine (Default)**: Rapid, zero-dependency analytical scoring based on the 11,365-guide benchmark weights for instant CPU screening (~150–600 ms per 1,000 guides on CPU, hardware-dependent; zero external model weight files required). Single-guide critiques (`explain_guide()`, `--explain`) evaluate biophysical rules derived from global TreeSHAP analysis.
> 2. **Full GBDT Engine (Optional)**: If trained or provided via `--model <path>` (e.g. generated via `python src/03_train_models.py --save_booster`), the predictor executes full multivariate gradient boosted tree inference across all 57 biophysical features.
> 3. **Global TreeSHAP Analysis**: Exact dataset-level Shapley value computation and 2-way interaction matrices are generated via `src/04_plot_shap_figures.py` and visualized in Figures 3A/3B.

---

## Repository Structure

```
.
├── src/
│   ├── 01_preprocess_data.py             # Data cleaning & Quality Control
│   ├── 01b_exploratory_data_analysis.py  # Exploratory Data Analysis & Figure 1 generator
│   ├── 01c_plot_figure_s1.py             # Supplementary Figure S1 generator
│   ├── 02_extract_features.py            # Master feature extraction pipeline
│   ├── 03_train_models.py                # 10-fold stratified CV & LightGBM training
│   ├── 04_plot_shap_figures.py           # TreeSHAP attribution & figure generator
│   ├── 05_model_validation_benchmark.py  # Benchmark validation & Figure 2
│   └── cas12a_xai/                       # Production Python package
│       ├── __init__.py
│       ├── features.py                   # SantaLucia thermodynamics & feature extraction
│       ├── predictor.py                  # Predictor class & biophysical critique engine
│       └── cli.py                        # Terminal CLI executable
├── tests/
│   ├── __init__.py
│   └── test_pipeline.py                  # Automated unit test suite
├── notebooks/
│   └── Cas12a_Complete_Master_Pipeline.ipynb # 1-Click Interactive Master Pipeline
├── data/
│   └── results/                          # Validated out-of-fold CV predictions & benchmark tables
├── requirements.txt
├── setup.py
├── pyproject.toml
├── LICENSE
├── CONTRIBUTING.md
└── README.md
```

---

## Running Automated Tests

Run the unit test suite across feature extractors, sequence validators, and predictors:
```bash
python -m unittest discover -s tests -p "test_*.py" -v
```

---

## Citation

If you use `cas12a-xai` or our biophysical findings in your research, please cite:

```bibtex
@article{tran2026cas12a_xai,
  title={Deciphering the Biophysical and Epistatic Determinants of CRISPR-Cas12a Guide RNA Cleavage Efficiency via Explainable Machine Learning},
  author={Tran, Minh},
  journal={Bioinformatics (In Review)},
  year={2026},
  url={https://github.com/minhtran-bio/cas12a-xai}
}
```

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.
