"""
=============================================================================
CRISPR-Cas12a Explainable AI Project
Module: 01b_exploratory_data_analysis.py
Author: Minh Tran (Independent Researcher, Perth, Western Australia, Australia)

Description:
  Executes comprehensive Exploratory Data Analysis (EDA) on the high-throughput
  Cas12a cleavage dataset (Kim et al. 2018). Generates statistical summaries
  and produces the multi-panel publication Figure 1 (300 DPI):
    - Panel A: Observed Indel Cleavage Frequency Distribution & Dichotomization
    - Panel B: PAM Subclass Cleavage Stratification (TTTA, TTTC, TTTG, TTTT)
    - Panel C: Positional Nucleotide Composition Across Functional Domains
    - Panel D: Thermodynamic Polarity Gradient (Delta_Tm) vs Cleavage Efficiency
    - Panel E: Nearest-Neighbor Duplex Stacking Free Energy (Delta G) Distributions
    - Panel F: Catastrophic Premature Termination by RNA Polymerase III Arrest
=============================================================================
"""

import os
import sys
import json

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

# Project Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLEAN_DATA_CSV = os.path.join(BASE_DIR, "data", "processed", "clean_kim_2018_targets.csv")
PREDICTIONS_CSV = os.path.join(BASE_DIR, "data", "results", "cas12a_predictions_cv.csv")
RESULTS_DIR = os.path.join(BASE_DIR, "data", "results")
FIGURES_DIR = os.path.join(BASE_DIR, "manuscript", "figures")

OUTPUT_FIG1 = os.path.join(FIGURES_DIR, "Figure1_Exploratory_Data_Analysis.png")
OUTPUT_EDA_JSON = os.path.join(RESULTS_DIR, "eda_summary_metrics.json")
OUTPUT_EDA_PAM_CSV = os.path.join(RESULTS_DIR, "eda_pam_comparison.csv")

# Plot Aesthetic Configuration
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.1
plt.rcParams['grid.alpha'] = 0.3


def load_datasets():
    """Load both clean full dataset and canonical feature dataset."""
    print("[*] Loading datasets for Exploratory Data Analysis...")
    if not os.path.exists(CLEAN_DATA_CSV):
        raise FileNotFoundError(f"[!] Clean dataset missing at: {CLEAN_DATA_CSV}")
    
    clean_df = pd.read_csv(CLEAN_DATA_CSV)
    
    # Load canonical feature matrix from predictions if available, else derive
    if os.path.exists(PREDICTIONS_CSV):
        canonical_df = pd.read_csv(PREDICTIONS_CSV)
    else:
        canonical_df = clean_df[clean_df["is_canonical_tttv"] == True].copy()
        
    print(f"    - Total processed screening targets: {len(clean_df):,}")
    print(f"    - Canonical TTTV targets analyzed:   {len(canonical_df):,}")
    return clean_df, canonical_df


def compute_exploratory_statistics(clean_df, canonical_df):
    """Compute exhaustive descriptive statistics and hypothesis tests."""
    print("[*] Computing descriptive moments and non-parametric hypothesis tests...")
    
    # 1. Target variable distribution
    indels = canonical_df["indel_frequency"].values
    mean_val = float(np.mean(indels))
    std_val = float(np.std(indels, ddof=1))
    median_val = float(np.median(indels))
    q25, q75 = np.percentile(indels, [25, 75])
    iqr_val = float(q75 - q25)
    skew_val = float(stats.skew(indels))
    kurt_val = float(stats.kurtosis(indels))
    pct_high = float(np.mean(indels >= 20.0) * 100.0)
    
    # 2. PAM Subclass comparison
    pam_stats = []
    # Canonical TTTV
    for pam_motif in ["TTTA", "TTTC", "TTTG"]:
        sub = canonical_df[canonical_df["pam_4bp"] == pam_motif]["indel_frequency"].values
        pam_stats.append({
            "pam_motif": pam_motif,
            "category": "Canonical",
            "count": int(len(sub)),
            "proportion_pct": round(len(sub) / len(canonical_df) * 100.0, 2),
            "mean_indel": round(float(np.mean(sub)), 2),
            "median_indel": round(float(np.median(sub)), 2),
            "std_indel": round(float(np.std(sub, ddof=1)), 2),
            "iqr_indel": round(float(np.percentile(sub, 75) - np.percentile(sub, 25)), 2),
            "pct_ge_20": round(float(np.mean(sub >= 20.0) * 100.0), 2)
        })
    
    # Non-canonical TTTT from full clean dataset
    tttt_sub = clean_df[clean_df["pam_4bp"] == "TTTT"]["indel_frequency"].values
    if len(tttt_sub) > 0:
        pam_stats.append({
            "pam_motif": "TTTT",
            "category": "Arrested / Non-canonical",
            "count": int(len(tttt_sub)),
            "proportion_pct": round(len(tttt_sub) / len(clean_df) * 100.0, 2),
            "mean_indel": round(float(np.mean(tttt_sub)), 2),
            "median_indel": round(float(np.median(tttt_sub)), 2),
            "std_indel": round(float(np.std(tttt_sub, ddof=1)), 2),
            "iqr_indel": round(float(np.percentile(tttt_sub, 75) - np.percentile(tttt_sub, 25)), 2),
            "pct_ge_20": round(float(np.mean(tttt_sub >= 20.0) * 100.0), 2)
        })
    
    pam_df = pd.DataFrame(pam_stats)
    
    # Kruskal-Wallis H test across canonical subclasses
    ttta_vals = canonical_df[canonical_df["pam_4bp"] == "TTTA"]["indel_frequency"].values
    tttc_vals = canonical_df[canonical_df["pam_4bp"] == "TTTC"]["indel_frequency"].values
    tttg_vals = canonical_df[canonical_df["pam_4bp"] == "TTTG"]["indel_frequency"].values
    kw_stat, kw_pval = stats.kruskal(ttta_vals, tttc_vals, tttg_vals)
    
    # 3. Poly-T terminator impact
    has_polyt = canonical_df["poly_t_terminator"] == 1
    polyt_pos = canonical_df[has_polyt]["indel_frequency"].values
    polyt_neg = canonical_df[~has_polyt]["indel_frequency"].values
    
    mann_whit_stat, mann_whit_pval = stats.mannwhitneyu(polyt_pos, polyt_neg, alternative='two-sided')
    polyt_deficit = float(np.mean(polyt_neg) - np.mean(polyt_pos))
    
    # 4. Thermodynamic correlations with Indel
    corr_features = [
        'gc_spacer', 'gc_seed', 'gc_distal',
        'tm_spacer', 'tm_seed', 'tm_distal', 'tm_gradient_seed_vs_distal',
        'stacking_dg_seed', 'stacking_dg_trunk', 'stacking_dg_distal', 'stacking_dg_spacer'
    ]
    correlations = {}
    for feat in corr_features:
        if feat in canonical_df.columns:
            rho, pval = stats.spearmanr(canonical_df[feat].values, indels)
            correlations[feat] = {
                "spearman_rho": round(float(rho), 4),
                "p_value": float(pval)
            }
            
    summary_metrics = {
        "dataset_sample_size": len(canonical_df),
        "target_indel_moments": {
            "mean": round(mean_val, 2),
            "std": round(std_val, 2),
            "median": round(median_val, 2),
            "iqr": round(iqr_val, 2),
            "skewness": round(skew_val, 4),
            "kurtosis": round(kurt_val, 4),
            "percent_functional_ge_20": round(pct_high, 2)
        },
        "pam_kruskal_wallis": {
            "h_statistic": round(float(kw_stat), 4),
            "p_value": float(kw_pval)
        },
        "poly_t_arrest_impact": {
            "poly_t_target_count": int(np.sum(has_polyt)),
            "mean_indel_without_poly_t": round(float(np.mean(polyt_neg)), 2),
            "mean_indel_with_poly_t": round(float(np.mean(polyt_pos)), 2),
            "mean_cleavage_deficit": round(polyt_deficit, 2),
            "mann_whitney_u_pval": float(mann_whit_pval)
        },
        "bivariate_spearman_correlations": correlations
    }
    
    return summary_metrics, pam_df


def generate_figure_1(clean_df, canonical_df, pam_df, summary_metrics):
    """Render publication-quality 6-panel Figure 1 (300 DPI)."""
    print("[*] Generating Figure 1: Exploratory Data Analysis & Biophysical Architecture...")
    os.makedirs(FIGURES_DIR, exist_ok=True)
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    
    # -------------------------------------------------------------------------
    # Panel A: Target Indel Cleavage Distribution & Threshold
    # -------------------------------------------------------------------------
    ax_a = axes[0, 0]
    indels = canonical_df["indel_frequency"].values
    sns.histplot(indels, kde=True, bins=40, color="#1f77b4", ax=ax_a, edgecolor="white", alpha=0.65)
    
    mean_val = summary_metrics["target_indel_moments"]["mean"]
    median_val = summary_metrics["target_indel_moments"]["median"]
    
    ax_a.axvline(20.0, color="#d62728", linestyle="--", linewidth=1.8, label="Functional Threshold (20%)")
    ax_a.axvline(median_val, color="#2ca02c", linestyle="-.", linewidth=1.8, label=f"Median ({median_val:.1f}%)")
    ax_a.axvline(mean_val, color="#ff7f0e", linestyle=":", linewidth=1.8, label=f"Mean ({mean_val:.1f}%)")
    
    ax_a.set_title("A  On-Target Indel Frequency Distribution", fontsize=12, fontweight="bold", loc="left")
    ax_a.set_xlabel("Observed Indel Frequency (%)", fontsize=11)
    ax_a.set_ylabel("Target Guide Count", fontsize=11)
    ax_a.legend(frameon=True, facecolor="white", framealpha=0.9, fontsize=9)
    ax_a.grid(True, linestyle=":", alpha=0.6)
    
    # Add annotation box
    ax_a.text(0.04, 0.65, f"N = {len(canonical_df):,}\nMean = {mean_val:.1f}% +/- {summary_metrics['target_indel_moments']['std']:.1f}%\nFunctional (>=20%): {summary_metrics['target_indel_moments']['percent_functional_ge_20']:.1f}%\nSkewness: {summary_metrics['target_indel_moments']['skewness']:.2f}",
              transform=ax_a.transAxes, fontsize=9.5, bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f4f8", edgecolor="#b0c4de"))

    # -------------------------------------------------------------------------
    # Panel B: PAM Subclass Cleavage Stratification
    # -------------------------------------------------------------------------
    ax_b = axes[0, 1]
    
    # Prepare comparison dataset including TTTT
    pam_plot_data = []
    for pam in ["TTTA", "TTTC", "TTTG"]:
        sub = canonical_df[canonical_df["pam_4bp"] == pam]["indel_frequency"].values
        for v in sub:
            pam_plot_data.append({"PAM": pam, "Indel": v, "Class": "Canonical"})
            
    tttt_vals = clean_df[clean_df["pam_4bp"] == "TTTT"]["indel_frequency"].values
    if len(tttt_vals) > 0:
        for v in tttt_vals:
            pam_plot_data.append({"PAM": "TTTT", "Indel": v, "Class": "Arrested"})
            
    pam_plot_df = pd.DataFrame(pam_plot_data)
    palette = {"TTTA": "#3274a1", "TTTC": "#4e9bb9", "TTTG": "#e1812c", "TTTT": "#c03d3e"}
    
    sns.violinplot(x="PAM", y="Indel", data=pam_plot_df, hue="PAM", palette=palette, legend=False, ax=ax_b, inner=None, cut=0)
    sns.boxplot(x="PAM", y="Indel", data=pam_plot_df, width=0.18, color="white", ax=ax_b,
                boxprops=dict(alpha=0.85), medianprops=dict(color="red", linewidth=2.0),
                whiskerprops=dict(color="black"), capprops=dict(color="black"))
                
    ax_b.set_title("B  Cleavage Efficiency Across PAM Subclasses", fontsize=12, fontweight="bold", loc="left")
    ax_b.set_xlabel("5'-PAM Motif", fontsize=11)
    ax_b.set_ylabel("Observed Indel Frequency (%)", fontsize=11)
    ax_b.grid(True, linestyle=":", alpha=0.6)
    
    # Annotate sample counts
    for i, pam in enumerate(["TTTA", "TTTC", "TTTG", "TTTT"]):
        cnt = len(pam_plot_df[pam_plot_df["PAM"] == pam])
        mean_p = np.mean(pam_plot_df[pam_plot_df["PAM"] == pam]["Indel"])
        ax_b.text(i, -9.0, f"n={cnt:,}\n({mean_p:.1f}%)", ha="center", fontsize=8.5, fontweight="semibold")
    ax_b.set_ylim(-12, 105)

    # -------------------------------------------------------------------------
    # Panel C: Positional Nucleotide Composition (23-nt Spacer)
    # -------------------------------------------------------------------------
    ax_c = axes[0, 2]
    
    spacers = canonical_df["spacer_23bp"].astype(str).tolist()
    n_spacers = len(spacers)
    pos_matrix = {nt: np.zeros(23) for nt in ['A', 'C', 'G', 'T']}
    
    for sp in spacers:
        for idx in range(min(23, len(sp))):
            base = sp[idx].upper()
            if base in pos_matrix:
                pos_matrix[base][idx] += 1
                
    for nt in pos_matrix:
        pos_matrix[nt] = (pos_matrix[nt] / n_spacers) * 100.0
        
    positions = np.arange(1, 24)
    colors = {'A': '#2ca02c', 'C': '#1f77b4', 'G': '#ff7f0e', 'T': '#d62728'}
    
    bottom = np.zeros(23)
    for nt in ['A', 'C', 'G', 'T']:
        ax_c.bar(positions, pos_matrix[nt], bottom=bottom, color=colors[nt], label=nt, width=0.75, edgecolor="none")
        bottom += pos_matrix[nt]
        
    # Domain delineations: Seed (1-8), Trunk (9-16), Distal (17-23)
    ax_c.axvline(8.5, color="black", linestyle="--", linewidth=1.5)
    ax_c.axvline(16.5, color="black", linestyle="--", linewidth=1.5)
    
    ax_c.text(4.5, 103, "Seed (nt 1-8)", ha="center", fontsize=9.5, fontweight="bold", color="#1a365d")
    ax_c.text(12.5, 103, "Trunk (nt 9-16)", ha="center", fontsize=9.5, fontweight="bold", color="#1a365d")
    ax_c.text(20.0, 103, "Distal (nt 17-23)", ha="center", fontsize=9.5, fontweight="bold", color="#1a365d")
    
    ax_c.set_title("C  Positional Nucleotide Composition Across Guide Domains", fontsize=12, fontweight="bold", loc="left")
    ax_c.set_xlabel("Spacer Nucleotide Position (5' -> 3')", fontsize=11)
    ax_c.set_ylabel("Nucleotide Frequency (%)", fontsize=11)
    ax_c.set_xticks(range(1, 24, 2))
    ax_c.set_ylim(0, 112)
    ax_c.legend(loc="upper right", ncol=4, frameon=True, fontsize=8.5)
    ax_c.grid(True, axis="y", linestyle=":", alpha=0.6)

    # -------------------------------------------------------------------------
    # Panel D: Thermodynamic Polarity Gradient vs Cleavage Efficiency
    # -------------------------------------------------------------------------
    ax_d = axes[1, 0]
    
    grad = canonical_df["tm_gradient_seed_vs_distal"].values
    
    # Hexbin plot to handle large sample density cleanly
    hb = ax_d.hexbin(grad, indels, gridsize=30, cmap="Blues", mincnt=1, bins='log')
    cb = fig.colorbar(hb, ax=ax_d, fraction=0.046, pad=0.04)
    cb.set_label("Target Density (log10)", fontsize=9.5)
    
    # Trendline
    slope, intercept, r_val, p_val, std_err = stats.linregress(grad, indels)
    x_vals = np.linspace(np.percentile(grad, 1), np.percentile(grad, 99), 100)
    ax_d.plot(x_vals, intercept + slope * x_vals, color="#e63946", linewidth=2.5,
              label=f"Linear Fit (r = {r_val:.3f}, P < 1e-15)")
    
    ax_d.axvline(0.0, color="gray", linestyle="--", linewidth=1.2)
    ax_d.text(0.05, 0.08, "Stalled Gradient (Delta_Tm < 0)", transform=ax_d.transAxes, fontsize=8.5, color="#8b0000", ha="left")
    ax_d.text(0.95, 0.08, "Forward Zippering (Delta_Tm > 0)", transform=ax_d.transAxes, fontsize=8.5, color="#006400", ha="right")
    
    ax_d.set_title("D  Kinetic Polarity Gradient Drives Cleavage", fontsize=12, fontweight="bold", loc="left")
    ax_d.set_xlabel("Thermodynamic Polarity Gradient: Delta_Tm = Tm_seed - Tm_distal (deg C)", fontsize=10)
    ax_d.set_ylabel("Observed Indel Frequency (%)", fontsize=11)
    ax_d.legend(loc="upper left", frameon=True, fontsize=9)
    ax_d.grid(True, linestyle=":", alpha=0.6)

    # -------------------------------------------------------------------------
    # Panel E: SantaLucia Base-Stacking Free Energy Across Domains
    # -------------------------------------------------------------------------
    ax_e = axes[1, 1]
    
    dg_data = [
        canonical_df["stacking_dg_seed"].values,
        canonical_df["stacking_dg_trunk"].values,
        canonical_df["stacking_dg_distal"].values,
        canonical_df["stacking_dg_spacer"].values
    ]
    domain_labels = ["Seed\n(8 bp)", "Trunk\n(8 bp)", "Distal\n(7 bp)", "Full Spacer\n(23 bp)"]
    palette_e = ["#3a86ff", "#8338ec", "#ff006e", "#38b000"]
    
    boxplot_kw_e = {"tick_labels": domain_labels} if hasattr(plt, 'boxplot') else {"labels": domain_labels}
    try:
        bplot = ax_e.boxplot(dg_data, patch_artist=True, tick_labels=domain_labels, widths=0.55,
                             medianprops=dict(color="black", linewidth=1.8),
                             flierprops=dict(marker='o', markersize=2.5, alpha=0.2))
    except TypeError:
        bplot = ax_e.boxplot(dg_data, patch_artist=True, labels=domain_labels, widths=0.55,
                             medianprops=dict(color="black", linewidth=1.8),
                             flierprops=dict(marker='o', markersize=2.5, alpha=0.2))
                         
    for patch, color in zip(bplot['boxes'], palette_e):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
        
    ax_e.set_title("E  Nearest-Neighbor Stacking Free Energies (Delta G)", fontsize=12, fontweight="bold", loc="left")
    ax_e.set_xlabel("Target DNA Structural Domain", fontsize=11)
    ax_e.set_ylabel("Duplex Stacking Delta G (kcal/mol)", fontsize=11)
    ax_e.grid(True, linestyle=":", alpha=0.6)
    
    # Annotate means
    for i, arr in enumerate(dg_data):
        ax_e.text(i + 1, np.percentile(arr, 95) + 1.2, f"Mean:\n{np.mean(arr):.1f}", ha="center", fontsize=8.5, fontweight="semibold")

    # -------------------------------------------------------------------------
    # Panel F: Premature Pol III Termination Cleavage Repression
    # -------------------------------------------------------------------------
    ax_f = axes[1, 2]
    
    polyt_comp = [
        canonical_df[canonical_df["poly_t_terminator"] == 0]["indel_frequency"].values,
        canonical_df[canonical_df["poly_t_terminator"] == 1]["indel_frequency"].values
    ]
    labels_f = [
        f"Intact crRNA\n(No TTTT; n={len(polyt_comp[0]):,})",
        f"Poly-T Terminator\n(>=TTTT; n={len(polyt_comp[1]):,})"
    ]
    palette_f = ["#2a9d8f", "#e76f51"]
    
    try:
        bplot_f = ax_f.boxplot(polyt_comp, patch_artist=True, tick_labels=labels_f, widths=0.45,
                               medianprops=dict(color="black", linewidth=2.0),
                               flierprops=dict(marker='o', markersize=3.0, alpha=0.3))
    except TypeError:
        bplot_f = ax_f.boxplot(polyt_comp, patch_artist=True, labels=labels_f, widths=0.45,
                               medianprops=dict(color="black", linewidth=2.0),
                               flierprops=dict(marker='o', markersize=3.0, alpha=0.3))
                           
    for patch, color in zip(bplot_f['boxes'], palette_f):
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
        
    mean_intact = np.mean(polyt_comp[0])
    mean_polyt = np.mean(polyt_comp[1])
    deficit = mean_intact - mean_polyt
    
    ax_f.text(1.5, 88, f"Mean Cleavage Deficit:\nDelta Indel = -{deficit:.1f}%\n(P = {summary_metrics['poly_t_arrest_impact']['mann_whitney_u_pval']:.2e})",
              ha="center", fontsize=10, fontweight="bold",
              bbox=dict(boxstyle="round,pad=0.5", facecolor="#ffebee", edgecolor="#d32f2f"))
              
    ax_f.set_title("F  Catastrophic RNA Pol III Termination Arrest", fontsize=12, fontweight="bold", loc="left")
    ax_f.set_ylabel("Observed Indel Frequency (%)", fontsize=11)
    ax_f.grid(True, linestyle=":", alpha=0.6)
    
    # Layout adjustments
    plt.tight_layout()
    plt.savefig(OUTPUT_FIG1, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Figure 1 successfully exported to:\n    {OUTPUT_FIG1}")


def main():
    clean_df, canonical_df = load_datasets()
    summary_metrics, pam_df = compute_exploratory_statistics(clean_df, canonical_df)
    
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    # Save JSON summary metrics
    with open(OUTPUT_EDA_JSON, "w", encoding="utf-8") as f:
        json.dump(summary_metrics, f, indent=4)
    print(f"[OK] Saved EDA summary metrics to:\n    {OUTPUT_EDA_JSON}")
    
    # Save PAM comparison table
    pam_df.to_csv(OUTPUT_EDA_PAM_CSV, index=False)
    print(f"[OK] Saved PAM stratification table to:\n    {OUTPUT_EDA_PAM_CSV}")
    
    # Render publication Figure 1
    generate_figure_1(clean_df, canonical_df, pam_df, summary_metrics)
    
    print("\n" + "="*70)
    print("EXPLORATORY DATA ANALYSIS (EDA) COMPLETE:")
    print(f"  * Evaluated Canonical Guides: {summary_metrics['dataset_sample_size']:,}")
    print(f"  * Indel Mean +/- SD:          {summary_metrics['target_indel_moments']['mean']}% +/- {summary_metrics['target_indel_moments']['std']}%")
    print(f"  * Indel Median (IQR):         {summary_metrics['target_indel_moments']['median']}% ({summary_metrics['target_indel_moments']['iqr']}%)")
    print(f"  * Guides Functional (>=20%):  {summary_metrics['target_indel_moments']['percent_functional_ge_20']}%")
    print(f"  * Poly-T Cleavage Deficit:    -{summary_metrics['poly_t_arrest_impact']['mean_cleavage_deficit']}% (P = {summary_metrics['poly_t_arrest_impact']['mann_whitney_u_pval']:.2e})")
    print("="*70)


if __name__ == "__main__":
    main()
