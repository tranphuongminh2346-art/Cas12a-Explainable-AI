"""
Feature Extraction and Biophysical Modeling for CRISPR-Cas12a Guide RNAs.

Implements SantaLucia (1998) unified nearest-neighbor thermodynamic parameters,
domain segmentation (Seed nt 1-8, Trunk nt 9-16, Distal nt 17-23),
directional kinetic polarity gradients, and positional seed encodings.
Total: Exactly 57 verified biophysical features.
"""

from typing import Dict, List, Union
import numpy as np
import pandas as pd

# SantaLucia 1998 Unified Nearest-Neighbor Base-Stacking Free Energies (kcal/mol at 37°C, 1M NaCl)
NN_DELTA_G: Dict[str, float] = {
    "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84
}

# Exactly 57 canonical features matching Supplementary Table S1
FEATURE_NAMES: List[str] = [
    'gc_34bp', 'gc_pam', 'gc_spacer', 'gc_seed', 'gc_trunk', 'gc_distal',
    'tm_spacer', 'tm_seed', 'tm_distal', 'tm_gradient_seed_vs_distal',
    'stacking_dg_seed', 'stacking_dg_trunk', 'stacking_dg_distal', 'stacking_dg_spacer',
    'poly_t_terminator', 'poly_t_strict_5t', 'poly_g_quadruplex',
    'pam_terminal_is_A', 'pam_terminal_is_C', 'pam_terminal_is_G', 'pam_terminal_is_T',
    'seed_dinuc_1_2_UU', 'seed_dinuc_1_2_CC', 'seed_dinuc_1_2_GG', 'seed_dinuc_1_2_GA'
] + [f'seed_pos_{p}_{b}' for p in range(1, 9) for b in ['A', 'C', 'G', 'T']]


def calc_gc(seq: str) -> float:
    """Calculate GC fraction of a DNA sequence."""
    if not seq:
        return 0.0
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq)


def calc_tm(seq: str) -> float:
    """
    Calculate melting temperature (°C) of a DNA/RNA duplex segment.
    Uses Wallace rule (<14 nt) or Marmur/Doty formula (>=14 nt).
    """
    if not seq:
        return 0.0
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    at = seq.count('A') + seq.count('T')
    if len(seq) < 14:
        return (at * 2.0) + (gc * 4.0)
    else:
        return 64.9 + 41.0 * (gc - 16.4) / len(seq)


def calc_stacking_dg(seq: str) -> float:
    """
    Calculate cumulative SantaLucia (1998) base stacking free energy (kcal/mol).
    Defaults unknown dinucleotide steps to mean Watson-Crick value (-1.20 kcal/mol).
    """
    if not seq or len(seq) < 2:
        return 0.0
    seq = seq.upper()
    return sum(NN_DELTA_G.get(seq[i:i+2], -1.20) for i in range(len(seq) - 1))


def decompose_target_sequence(sequence: str) -> Dict[str, str]:
    """
    Parse an input sequence into coordinate domains.
    Accepts:
      - 34-bp target site (4bp 5' flank + 4bp PAM + 20bp guide + 3bp target context + 3bp 3' flank)
      - 27-bp PAM + 23-bp protospacer window (4bp PAM + 20bp guide + 3bp target context)
      - 24-bp PAM + 20-nt guide (4bp PAM + 20bp guide; target context neutrally padded)
      - 23-bp protospacer window (20bp guide + 3bp target context; assumes TTTC PAM)
      - 20-bp guide alone (assumes canonical TTTC PAM and neutral context padding)
    """
    s = sequence.strip().upper()
    if len(s) == 34:
        upstream = s[0:4]
        pam = s[4:8]
        spacer = s[8:31]
        downstream = s[31:34]
    elif len(s) == 27:
        upstream = "NNNN"
        pam = s[0:4]
        spacer = s[4:27]
        downstream = "NNN"
    elif len(s) == 24:
        upstream = "NNNN"
        pam = s[0:4]
        spacer = s[4:24] + "NNN"
        downstream = "NNN"
    elif len(s) == 23:
        upstream = "NNNN"
        pam = "TTTC"  # Default canonical
        spacer = s
        downstream = "NNN"
    elif len(s) == 20:
        upstream = "NNNN"
        pam = "TTTC"  # Default canonical
        spacer = s + "NNN"
        downstream = "NNN"
    else:
        raise ValueError(
            f"Invalid sequence length {len(s)} bp for '{s}'. Expected 34 bp, 27 bp, 24 bp, 23 bp, or 20 bp."
        )

    seed = spacer[0:8]
    trunk = spacer[8:16]
    distal = spacer[16:23]
    guide_20bp = spacer[0:20]
    target_context_3bp = spacer[20:23]

    has_guide20_polyt = int("TTTT" in guide_20bp)
    has_boundary_polyt = int(("TTTT" in spacer) and ("TTTT" not in guide_20bp))

    return {
        "raw_sequence": s,
        "upstream_4bp": upstream,
        "pam_4bp": pam,
        "spacer_23bp": spacer,
        "guide_20bp": guide_20bp,
        "target_context_3bp": target_context_3bp,
        "has_guide20_polyt": has_guide20_polyt,
        "has_boundary_polyt": has_boundary_polyt,
        "seed_8bp": seed,
        "trunk_8bp": trunk,
        "distal_7bp": distal,
        "downstream_3bp": downstream,
    }


def extract_features_single(sequence: str) -> Dict[str, Union[float, int]]:
    """
    Extract full 57-dimensional biophysical feature vector from a single target sequence.
    """
    d = decompose_target_sequence(sequence)
    s34 = d["raw_sequence"] if len(d["raw_sequence"]) == 34 else f"{d['upstream_4bp']}{d['pam_4bp']}{d['spacer_23bp']}{d['downstream_3bp']}"
    pam = d["pam_4bp"]
    sp = d["spacer_23bp"]
    seed = d["seed_8bp"]
    trunk = d["trunk_8bp"]
    distal = d["distal_7bp"]

    tm_s = calc_tm(seed)
    tm_d = calc_tm(distal)
    d12 = seed[0:2] if len(seed) >= 2 else ""

    features: Dict[str, Union[float, int]] = {
        'gc_34bp': calc_gc(s34),
        'gc_pam': calc_gc(pam),
        'gc_spacer': calc_gc(sp),
        'gc_seed': calc_gc(seed),
        'gc_trunk': calc_gc(trunk),
        'gc_distal': calc_gc(distal),
        'tm_spacer': calc_tm(sp),
        'tm_seed': tm_s,
        'tm_distal': tm_d,
        'tm_gradient_seed_vs_distal': tm_s - tm_d,
        'stacking_dg_seed': calc_stacking_dg(seed),
        'stacking_dg_trunk': calc_stacking_dg(trunk),
        'stacking_dg_distal': calc_stacking_dg(distal),
        'stacking_dg_spacer': calc_stacking_dg(sp),
        'poly_t_terminator': int('TTTT' in sp),
        'poly_t_strict_5t': int('TTTTT' in sp),
        'poly_g_quadruplex': int('GGGG' in sp),
        'pam_terminal_is_A': int(pam[-1:] == 'A'),
        'pam_terminal_is_C': int(pam[-1:] == 'C'),
        'pam_terminal_is_G': int(pam[-1:] == 'G'),
        'pam_terminal_is_T': int(pam[-1:] == 'T'),
        'seed_dinuc_1_2_UU': int(d12 == 'TT'),
        'seed_dinuc_1_2_CC': int(d12 == 'CC'),
        'seed_dinuc_1_2_GG': int(d12 == 'GG'),
        'seed_dinuc_1_2_GA': int(d12 == 'GA'),
    }

    # One-hot encoding for seed positions 1 to 8 (32 features)
    for p in range(8):
        base_char = seed[p] if p < len(seed) else ''
        for base in ['A', 'C', 'G', 'T']:
            features[f'seed_pos_{p+1}_{base}'] = int(base_char == base)

    return features


def extract_features_batch(sequences: List[str]) -> pd.DataFrame:
    """
    Extract biophysical feature matrix for a collection of sequences.
    Returns a pandas DataFrame with columns strictly ordered according to FEATURE_NAMES (57 features).
    """
    records = [extract_features_single(seq) for seq in sequences]
    df = pd.DataFrame(records)
    return df[FEATURE_NAMES]
