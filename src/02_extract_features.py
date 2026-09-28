"""
=============================================================================
CRISPR-Cas12a Explainable AI Project
Module: 02_extract_features.py
Author: Minh Tran (School of Molecular Sciences, University of Western Australia)

Description:
  Computes comprehensive biophysical, sequence, and thermodynamic features
  from clean Cas12a target sequences:
    1. Regional GC contents (PAM, Seed, Trunk, Distal, Flanks)
    2. Nearest-Neighbor Thermodynamics (Melting Temp Tm, Stacking Energy dG)
    3. Thermodynamic Polarity Gradient (Seed Tm vs Distal Tm)
    4. Sequence Motif Penalties (Poly-T Pol III terminators, G-quadruplexes)
    5. Positional k-mer encodings across the 23bp protospacer and PAM
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

import pandas as pd
import numpy as np

# Project Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INPUT_CLEAN_CSV = os.path.join(BASE_DIR, "data", "processed", "clean_kim_2018_targets.csv")
OUTPUT_FEATURE_CSV = os.path.join(BASE_DIR, "data", "processed", "cas12a_features_master.csv")

# SantaLucia 1998 / Sugimoto 1996 unified nearest-neighbor thermodynamic parameters
# Delta G (kcal/mol at 37°C, 1M NaCl) for dinucleotide stacking
NN_DELTA_G = {
    "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84
}

def calc_gc(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)

def calc_tm(seq: str) -> float:
    """Nearest-neighbor / Wallace empirical melting temperature."""
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    if len(seq) < 14:
        return (at * 2.0) + (gc * 4.0)
    else:
        return 64.9 + 41.0 * (gc - 16.4) / len(seq)

def calc_stacking_energy(seq: str) -> float:
    """Cumulative nearest-neighbor stacking free energy (kcal/mol). More negative = more stable."""
    if len(seq) < 2:
        return 0.0
    total_dg = 0.0
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        total_dg += NN_DELTA_G.get(dinuc, -1.20)
    return total_dg

def extract_features_for_row(row):
    seq34 = str(row["raw_sequence"]).upper()
    upstream = seq34[0:4]
    pam = seq34[4:8]
    spacer = seq34[8:31]
    downstream = seq34[31:34]
    
    # 4 Functional Kinetic Domains (Strohkendl et al. Molecular Cell 2018)
    seed = spacer[0:8]       # Positions 1 to 8: R-loop initiation
    trunk = spacer[8:16]     # Positions 9 to 16: Conformational locking
    distal = spacer[16:23]   # Positions 17 to 23: Cleavage kinetics fine-tuning
    
    tm_spacer = calc_tm(spacer)
    tm_seed = calc_tm(seed)
    tm_distal = calc_tm(distal)
    
    feat = {
        # 1. GC Contents across functional domains
        "gc_34bp": calc_gc(seq34),
        "gc_upstream": calc_gc(upstream),
        "gc_pam": calc_gc(pam),
        "gc_spacer": calc_gc(spacer),
        "gc_seed": calc_gc(seed),
        "gc_trunk": calc_gc(trunk),
        "gc_distal": calc_gc(distal),
        "gc_downstream": calc_gc(downstream),
        
        # 2. Thermodynamics & Stability
        "tm_spacer": tm_spacer,
        "tm_seed": tm_seed,
        "tm_distal": tm_distal,
        "tm_gradient_seed_vs_distal": tm_seed - tm_distal, # Thermodynamic polarity
        "stacking_dg_seed": calc_stacking_energy(seed),
        "stacking_dg_trunk": calc_stacking_energy(trunk),
        "stacking_dg_distal": calc_stacking_energy(distal),
        "stacking_dg_spacer": calc_stacking_energy(spacer),
        
        # 3. Biological Sequence Penalties (DeWeirdt et al. 2021)
        "poly_t_terminator": int("TTTT" in spacer),
        "poly_t_strict_5t": int("TTTTT" in spacer),
        "poly_g_quadruplex": int("GGGG" in spacer),
        "poly_c_count": int("CCCC" in spacer),
        
        # 4. Terminal PAM base (Position -1)
        "pam_terminal_is_A": int(pam[3] == "A"),
        "pam_terminal_is_C": int(pam[3] == "C"),
        "pam_terminal_is_G": int(pam[3] == "G"),
        "pam_terminal_is_T": int(pam[3] == "T") # Disfavored in canonical AsCas12a
    }
    
    # 5. Positional one-hot encoding for seed region (positions 1-8)
    for pos, base in enumerate(seed, start=1):
        for b in ["A", "C", "G", "T"]:
            feat[f"seed_pos_{pos}_{b}"] = int(base == b)
            
    # 6. Selected key dinucleotides in seed (positions 1-2, 2-3)
    feat["seed_dinuc_1_2_UU"] = int(seed[0:2] == "TT") # Strongly disfavored in CRISPR-DT
    feat["seed_dinuc_1_2_CC"] = int(seed[0:2] == "CC")
    feat["seed_dinuc_1_2_GG"] = int(seed[0:2] == "GG")
    feat["seed_dinuc_1_2_GA"] = int(seed[0:2] == "GA")
    
    return feat

def main():
    print(f"[*] Reading clean targets from: {INPUT_CLEAN_CSV}")
    if not os.path.exists(INPUT_CLEAN_CSV):
        raise FileNotFoundError(f"[!] Please run 01_preprocess_data.py first!")
        
    clean_df = pd.read_csv(INPUT_CLEAN_CSV)
    print(f"[+] Loaded {len(clean_df)} target sequences.")
    
    print("[*] Computing biophysical and thermodynamic features...")
    feature_list = []
    for idx, row in clean_df.iterrows():
        feat = extract_features_for_row(row)
        feature_list.append(feat)
        if (idx + 1) % 4000 == 0 or (idx + 1) == len(clean_df):
            print(f"    - Processed {idx + 1}/{len(clean_df)} sequences...")
            
    feat_df = pd.DataFrame(feature_list)
    
    # Concatenate targets metadata with features
    final_df = pd.concat([
        clean_df[["raw_sequence", "pam_4bp", "spacer_23bp", "indel_frequency", "is_canonical_tttv", "is_high_efficiency"]],
        feat_df
    ], axis=1)
    
    print(f"[OK] Feature extraction complete. Total feature columns: {feat_df.shape[1]}")
    final_df.to_csv(OUTPUT_FEATURE_CSV, index=False)
    print(f"[OK] Master feature matrix saved to:\n    {OUTPUT_FEATURE_CSV}")

if __name__ == "__main__":
    main()
