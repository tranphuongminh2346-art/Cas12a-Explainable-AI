"""
Feature Extraction Pipeline for CRISPR-Cas12a Guide RNA Sequences
Extracts:
  1. Positional one-hot and k-mer encodings (1-mer, 2-mer).
  2. Sub-region GC contents (PAM, Seed, Trunk, Distal).
  3. Homopolymer counts (e.g. poly-T terminator check).
  4. Thermodynamic and hybrid stability indicators (melting temperature Tm, delta G approximations).
"""

import numpy as np
import pandas as pd

# Nearest-neighbor thermodynamic parameters (SantaLucia 1998 / Sugimoto 1995 approximation)
# Delta G (kcal/mol) for dinucleotide pairs
NN_DELTA_G = {
    "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84
}

def calculate_gc(seq: str) -> float:
    if not seq:
        return 0.0
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)

def approximate_tm(seq: str) -> float:
    """Basic nearest-neighbor / Wallace rule melting temperature approximation."""
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    if len(seq) < 14:
        return (at * 2) + (gc * 4)
    else:
        return 64.9 + 41.0 * (gc - 16.4) / len(seq)

def calculate_stacking_energy(seq: str) -> float:
    """Compute cumulative nearest-neighbor stacking free energy (kcal/mol)."""
    seq = seq.upper()
    total_dg = 0.0
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        total_dg += NN_DELTA_G.get(dinuc, -1.20) # default fallback
    return total_dg

def extract_cas12a_features(target_34bp: str) -> dict:
    """
    Extract comprehensive biophysical & sequence features from a 34bp target:
      - 4bp upstream: positions 0..3
      - 4bp PAM (TTTV): positions 4..7
      - 23bp protospacer: positions 8..30
        - Seed (pos 8..15, length 8)
        - Mid-trunk (pos 16..23, length 8)
        - Distal (pos 24..30, length 7)
      - 3bp downstream: positions 31..33
    """
    seq = target_34bp.upper().replace("U", "T")
    assert len(seq) == 34, f"Target sequence must be exactly 34 bp, got {len(seq)}"
    
    upstream = seq[0:4]
    pam = seq[4:8]
    protospacer = seq[8:31]
    downstream = seq[31:34]
    
    seed = protospacer[0:8]
    trunk = protospacer[8:16]
    distal = protospacer[16:23]
    
    features = {
        # Overall regional GC contents
        "gc_total_34bp": calculate_gc(seq),
        "gc_pam": calculate_gc(pam),
        "gc_protospacer": calculate_gc(protospacer),
        "gc_seed": calculate_gc(seed),
        "gc_trunk": calculate_gc(trunk),
        "gc_distal": calculate_gc(distal),
        
        # Thermodynamics & hybrid stability
        "tm_protospacer": approximate_tm(protospacer),
        "tm_seed": approximate_tm(seed),
        "tm_distal": approximate_tm(distal),
        "stacking_dg_seed": calculate_stacking_energy(seed),
        "stacking_dg_protospacer": calculate_stacking_energy(protospacer),
        
        # Sequence penalties & structural signals
        "poly_t_count": int("TTTT" in protospacer),
        "poly_g_count": int("GGGG" in protospacer),
        "pam_terminal_base": pam[-1], # Should be A, C, or G (V base)
        "pam_is_valid_tttv": int(pam.startswith("TTT") and pam[3] in ["A", "C", "G"])
    }
    
    # Positional nucleotides in seed region (positions 1-8 of protospacer)
    for idx, char in enumerate(seed, start=1):
        for base in ["A", "C", "G", "T"]:
            features[f"pos_seed_{idx}_{base}"] = int(char == base)
            
    return features

if __name__ == "__main__":
    test_seq = "ACCGTTTAAGCTAGCTAGCTAGCTAGCTAGCTAAG"
    feats = extract_cas12a_features(test_seq)
    print(f"[OK] Feature extraction test successful. Total features generated: {len(feats)}")
