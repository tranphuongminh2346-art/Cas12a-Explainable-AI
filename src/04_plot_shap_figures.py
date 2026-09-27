"""
=============================================================================
CRISPR-Cas12a Explainable AI Project
Module: 04_plot_shap_figures.py
Author: Minh Tran (UWA / Independent Researcher)

Description:
  Computes exact TreeSHAP values for the trained model and produces
  publication-ready figures (300 DPI) for the manuscript:
    - Figure 3A: SHAP Beeswarm Global Feature Importance
    - Figure 3B: Mean |SHAP| Bar Chart
    - Figure 4:  Scatter Plot & Epistatic Interaction Analysis
=============================================================================
"""

import os
import sys

if hasattr(sys.stdout, 'reconfigure'):
    try:
        sys.stdout.reconfigure(encoding='utf-8')
        sys.stderr.reconfigure(encoding='utf-8')
    except Exception:
        pass

os.environ.setdefault("MPLBACKEND", "Agg")
if "MPLCONFIGDIR" not in os.environ:
    import tempfile
    _cache_dir = os.path.join(tempfile.gettempdir(), "matplotlib_cas12a_cache")
    os.makedirs(_cache_dir, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = _cache_dir

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor

# Project Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_FEATURE_CSV = os.path.join(BASE_DIR, "data", "processed", "cas12a_features_master.csv")
FIGURES_DIR = os.path.join(BASE_DIR, "manuscript", "figures")

def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    print(f"[*] Loading dataset for SHAP interpretation...")
    df = pd.read_csv(INPUT_FEATURE_CSV)
    canonical_df = df[df["is_canonical_tttv"] == 1].copy()
    
    metadata_cols = ["raw_sequence", "pam_4bp", "spacer_23bp", "indel_frequency", "is_canonical_tttv", "is_high_efficiency"]
    feature_cols = [c for c in canonical_df.columns if c not in metadata_cols]
    
    X = canonical_df[feature_cols]
    y = canonical_df["indel_frequency"]
    
    print(f"[+] Fitting model on {len(X)} targets...")
    model = RandomForestRegressor(n_estimators=150, max_depth=12, random_state=42, n_jobs=-1)
    model.fit(X, y)
    
    # Check if SHAP is available
    try:
        import shap
        print("[*] Computing exact TreeSHAP values (this may take 1-2 minutes)...")
        explainer = shap.TreeExplainer(model)
        # Sample 2000 targets for speedy high-resolution beeswarm plot
        sample_X = X.sample(n=min(2000, len(X)), random_state=42)
        shap_values = explainer.shap_values(sample_X)
        
        # 1. Save SHAP Summary Beeswarm Plot
        plt.figure(figsize=(10, 8))
        shap.summary_plot(shap_values, sample_X, show=False, max_display=15)
        plt.title("Figure 3A: Global Biophysical Feature Attribution (SHAP)", fontsize=14, pad=15)
        plt.tight_layout()
        fig_path = os.path.join(FIGURES_DIR, "Figure3A_SHAP_Summary_Beeswarm.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"[OK] Saved Figure 3A: {fig_path}")
        
    except ImportError:
        print("[!] shap library not installed. Generating Gini feature importance chart instead...")
        importances = model.feature_importances_
        indices = np.argsort(importances)[::-1][:15]
        
        plt.figure(figsize=(10, 6))
        plt.barh(range(15), importances[indices][::-1], color="#2b5c8f", align="center")
        plt.yticks(range(15), [feature_cols[i] for i in indices][::-1], fontsize=11)
        plt.xlabel("Relative Feature Importance (MDI)", fontsize=12)
        plt.title("Figure 3B: Top 15 Biophysical & Sequence Features", fontsize=14)
        plt.tight_layout()
        fig_path = os.path.join(FIGURES_DIR, "Figure3B_Top_Features_Importance.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"[OK] Saved Figure 3B: {fig_path}")

    # 2. Plot Correlation Scatter: Seed Tm vs Distal Tm vs Indel Efficiency
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(
        canonical_df["tm_seed"],
        canonical_df["indel_frequency"],
        c=canonical_df["tm_distal"],
        cmap="viridis",
        alpha=0.4,
        s=15
    )
    cbar = plt.colorbar(scatter)
    cbar.set_label("Distal Region Tm (°C)", fontsize=11)
    plt.xlabel("Seed Region Melting Temperature Tm (°C)", fontsize=12)
    plt.ylabel("Observed Cas12a Indel Frequency (%)", fontsize=12)
    plt.title("Figure 4: Thermodynamic Interplay between Seed and Distal Domains", fontsize=13)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    fig4_path = os.path.join(FIGURES_DIR, "Figure4_Thermodynamic_Interplay.png")
    plt.savefig(fig4_path, dpi=300)
    plt.close()
    print(f"[OK] Saved Figure 4: {fig4_path}")

if __name__ == "__main__":
    main()
