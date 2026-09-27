"""
=============================================================================
CRISPR-Cas12a Explainable AI Project
Module: 01c_plot_figure_s1.py
Author: Minh Tran (UWA / Independent Researcher)

Description:
  Generates publication-quality Supplementary Figure S1 (300 DPI) adhering
  strictly to Bioinformatics / Nature Biotechnology visual formatting:
    - Panel A: Positional Nucleotide Enrichment (Log2 Odds Ratio) across 23-nt Protospacer
               with dedicated dual-axis domain annotation (Seed, Trunk, Distal)
    - Panel B: Lower-Triangular Bivariate Spearman Rank Correlation Heatmap
    - Panel C: Cleavage Efficiency Distribution across Stacking Free Energy (dG)
               Quartiles for Seed, Trunk, and Distal Domains
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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PREDICTIONS_CSV = os.path.join(BASE_DIR, "data", "results", "cas12a_predictions_cv.csv")
FIGURES_DIR = os.path.join(BASE_DIR, "manuscript", "figures")
PACKAGE_FIG_DIR = os.path.join(BASE_DIR, "cas12a_10fold_results_package", "figures")

OUTPUT_FIG_S1 = os.path.join(FIGURES_DIR, "Figure_S1_Biophysical_Correlation_Matrix.png")
OUTPUT_FIG_S1_PKG = os.path.join(PACKAGE_FIG_DIR, "Figure_S1_Biophysical_Correlation_Matrix.png")

# Aesthetic Configuration matching Figure 2 and Figure S2
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['axes.edgecolor'] = '#1a1a1a'
plt.rcParams['grid.alpha'] = 0.4
plt.rcParams['grid.linestyle'] = ':'


def compute_log2_odds_ratio(df):
    """Compute positional nucleotide log2 odds ratio for >=20% vs <20% indels."""
    high = df[df["indel_frequency"] >= 20.0]["spacer_23bp"].tolist()
    low = df[df["indel_frequency"] < 20.0]["spacer_23bp"].tolist()
    
    n_high = len(high)
    n_low = len(low)
    
    nucleotides = ['A', 'C', 'G', 'T']
    odds_ratio_matrix = np.zeros((4, 23))
    
    for pos in range(23):
        for idx, nt in enumerate(nucleotides):
            cnt_high = sum(1 for s in high if len(s) > pos and s[pos].upper() == nt)
            cnt_low = sum(1 for s in low if len(s) > pos and s[pos].upper() == nt)
            
            # Laplace smoothing (pseudocount = 1)
            p_high = (cnt_high + 1) / (n_high + 4)
            p_low = (cnt_low + 1) / (n_low + 4)
            
            odds_high = p_high / (1.0 - p_high)
            odds_low = p_low / (1.0 - p_low)
            
            odds_ratio_matrix[idx, pos] = np.log2(odds_high / odds_low)
            
    return odds_ratio_matrix, nucleotides


def generate_figure_s1():
    print("[*] Generating Supplementary Figure S1: Protospacer Odds Ratios & Biophysical Correlation Atlas...")
    if not os.path.exists(PREDICTIONS_CSV):
        raise FileNotFoundError(f"[!] Predictions CSV not found at: {PREDICTIONS_CSV}")
        
    df = pd.read_csv(PREDICTIONS_CSV)
    print(f"    - Loaded N = {len(df):,} canonical Cas12a target guides.")
    
    os.makedirs(FIGURES_DIR, exist_ok=True)
    os.makedirs(PACKAGE_FIG_DIR, exist_ok=True)
    
    fig = plt.figure(figsize=(19, 13))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.85, 1.15], hspace=0.38, wspace=0.25)
    
    # =========================================================================
    # Panel A: Positional Nucleotide Log2 Odds Ratio Heatmap (Spans full top row)
    # =========================================================================
    ax_a = fig.add_subplot(gs[0, :])
    
    odds_matrix, nucleotides = compute_log2_odds_ratio(df)
    vbound = max(abs(odds_matrix.min()), abs(odds_matrix.max()))
    vbound = round(float(vbound) + 0.05, 1)
    
    sns.heatmap(
        odds_matrix,
        cmap="coolwarm",
        vmin=-vbound,
        vmax=vbound,
        annot=True,
        fmt=".2f",
        annot_kws={"size": 9.0, "weight": "normal"},
        cbar_kws={"label": "Log2 Odds Ratio (Indel >= 20% vs < 20%)", "shrink": 0.90, "pad": 0.02},
        yticklabels=nucleotides,
        xticklabels=[str(i) for i in range(1, 24)],
        ax=ax_a,
        linewidths=0.6,
        linecolor="#ffffff"
    )
    
    # Domain boundary vertical lines
    ax_a.axvline(8, color="#111827", linestyle="--", linewidth=1.8, alpha=0.9)
    ax_a.axvline(16, color="#111827", linestyle="--", linewidth=1.8, alpha=0.9)
    
    # Secondary top X-axis for domain annotations: avoids ANY vertical collision with title
    ax_top = ax_a.twiny()
    ax_top.set_xlim(ax_a.get_xlim())
    ax_top.set_xticks([4.0, 12.0, 19.5])
    ax_top.set_xticklabels([
        "Seed Domain (Positions 1–8)",
        "Trunk Domain (Positions 9–16)",
        "Distal Domain (Positions 17–23)"
    ], fontsize=11, fontweight="bold", color="#1e3a8a")
    ax_top.tick_params(axis='x', length=0, pad=8)
    
    # Elevated main title cleanly above secondary axis
    ax_a.set_title("A. Positional Nucleotide Enrichment Across Protospacer (Log2 Odds Ratio: Functional >= 20% vs < 20%)",
                   fontsize=13, fontweight="bold", loc="left", pad=32)
    ax_a.set_xlabel("Protospacer Nucleotide Position (5' to 3')", fontsize=11, fontweight="bold", labelpad=8)
    ax_a.set_ylabel("Base", fontsize=11, fontweight="bold")
    ax_a.tick_params(axis='x', labelsize=10)
    ax_a.tick_params(axis='y', labelsize=11, rotation=0)
    
    # =========================================================================
    # Panel B: Spearman Bivariate Correlation Matrix
    # =========================================================================
    ax_b = fig.add_subplot(gs[1, 0])
    
    corr_features = [
        ('indel_frequency', 'Indel (%)'),
        ('gc_spacer', 'Spacer GC'),
        ('gc_seed', 'Seed GC'),
        ('gc_distal', 'Distal GC'),
        ('tm_seed', 'Seed Tm'),
        ('tm_distal', 'Distal Tm'),
        ('tm_gradient_seed_vs_distal', 'Polarity Delta Tm'),
        ('stacking_dg_spacer', 'Spacer dG'),
        ('stacking_dg_seed', 'Seed dG'),
        ('stacking_dg_distal', 'Distal dG')
    ]
    
    col_names = [f[0] for f in corr_features]
    display_names = [f[1] for f in corr_features]
    
    sub_df = df[col_names].copy()
    corr_matrix = sub_df.corr(method="spearman").values
    
    # Mask strictly upper triangle
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    
    sns.heatmap(
        corr_matrix,
        mask=mask,
        cmap="vlag",
        vmin=-1.0,
        vmax=1.0,
        annot=True,
        fmt=".2f",
        annot_kws={"size": 8.5},
        cbar_kws={"label": "Spearman Rank Correlation (rho)", "shrink": 0.85, "pad": 0.03},
        xticklabels=display_names,
        yticklabels=display_names,
        ax=ax_b,
        linewidths=0.6,
        linecolor="#ffffff"
    )
    
    ax_b.set_title("B. Bivariate Spearman Correlation Matrix of Biophysical Parameters",
                   fontsize=12, fontweight="bold", loc="left", pad=12)
    ax_b.tick_params(axis='x', rotation=45, labelsize=9.5)
    ax_b.tick_params(axis='y', rotation=0, labelsize=9.5)
    
    # =========================================================================
    # Panel C: Stacking Free Energy Quartiles vs Indel Efficiency
    # =========================================================================
    ax_c = fig.add_subplot(gs[1, 1])
    
    # Stratify Seed, Trunk, and Distal stacking energies into 4 quartiles
    # In SantaLucia stacking, more negative dG = higher stability
    plot_rows = []
    domains = [
        ("Seed (nt 1–8)", "stacking_dg_seed"),
        ("Trunk (nt 9–16)", "stacking_dg_trunk"),
        ("Distal (nt 17–23)", "stacking_dg_distal")
    ]
    
    for dom_label, col in domains:
        vals = df[col].values
        # Negative values: lower value means more negative = higher stability
        # We sort quartiles such that Q4 = Most negative (highest stability), Q1 = Least negative (lowest stability)
        q_bins = pd.qcut(vals, q=4, labels=["Q4 (Max Stability)", "Q3", "Q2", "Q1 (Min Stability)"])
        q_bins_str = pd.Series(q_bins).astype(str)

        for q, indel in zip(q_bins_str, df["indel_frequency"].values):
            plot_rows.append({
                "Domain": dom_label,
                "Quartile": q,
                "Indel": indel
            })
            
    q_df = pd.DataFrame(plot_rows)
    order_q = ["Q1 (Min Stability)", "Q2", "Q3", "Q4 (Max Stability)"]
    palette_c = ["#93c5fd", "#3b82f6", "#1d4ed8", "#1e3a8a"]
    
    sns.boxplot(
        x="Domain",
        y="Indel",
        hue="Quartile",
        hue_order=order_q,
        data=q_df,
        palette=palette_c,
        ax=ax_c,
        width=0.72,
        fliersize=1.5,
        flierprops=dict(alpha=0.15, marker='o'),
        medianprops=dict(color="#dc2626", linewidth=1.8),
        boxprops=dict(alpha=0.88)
    )
    
    ax_c.axhline(20.0, color="#b91c1c", linestyle="--", linewidth=1.4, alpha=0.85, label="Functional Threshold (20%)")
    ax_c.axhline(52.1, color="#4b5563", linestyle=":", linewidth=1.4, alpha=0.85, label="Cohort Mean Indel (52.1%)")
    
    ax_c.set_title("C. Cleavage Efficiency Across Stacking Stability Quartiles by Domain",
                   fontsize=12, fontweight="bold", loc="left", pad=12)
    ax_c.set_xlabel("Target DNA Functional Domain", fontsize=11, fontweight="bold", labelpad=8)
    ax_c.set_ylabel("Observed Cas12a Indel Frequency (%)", fontsize=11, fontweight="bold", labelpad=8)
    ax_c.set_ylim(-5, 108)
    ax_c.grid(True, axis='y', linestyle=':', alpha=0.6)
    ax_c.tick_params(axis='both', labelsize=10)
    
    # Legend formatting
    handles, labels = ax_c.get_legend_handles_labels()
    ax_c.legend(handles=handles, labels=labels, loc="lower left", frameon=True,
                facecolor="#ffffff", edgecolor="#d1d5db", fontsize=8.8, ncol=2)
    
    # Clean professional margins without tight_layout collision
    plt.subplots_adjust(top=0.92, bottom=0.08, left=0.06, right=0.95, hspace=0.38, wspace=0.25)
    
    plt.savefig(OUTPUT_FIG_S1, dpi=300)
    plt.savefig(OUTPUT_FIG_S1_PKG, dpi=300)
    plt.close()
    
    print(f"[OK] Saved Supplementary Figure S1 to: {OUTPUT_FIG_S1}")
    print(f"[OK] Mirrored to: {OUTPUT_FIG_S1_PKG}")


if __name__ == "__main__":
    generate_figure_s1()
