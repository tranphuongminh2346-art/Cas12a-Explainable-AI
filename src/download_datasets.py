"""
Dataset Acquisition Script for CRISPR-Cas12a Benchmark Data
Datasets:
  1. Kim et al. (2018) Nature Biotechnology (DeepCpf1 HT1, HT2, HT3, HEK-lenti, HEK-plasmid)
  2. DeWeirdt et al. (2021) Nature Biotechnology (AsCas12a combinatorial screen)
"""

import os
import urllib.request
import pandas as pd

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
RAW_DIR = os.path.join(DATA_DIR, "raw")
PROCESSED_DIR = os.path.join(DATA_DIR, "processed")

def setup_directories():
    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(PROCESSED_DIR, exist_ok=True)
    print(f"[+] Directories verified:\n    Raw: {RAW_DIR}\n    Processed: {PROCESSED_DIR}")

def fetch_kim_supplementary():
    """
    Kim et al. 2018 Supplementary Data contains the ~16,292 HT1 targets with indel frequencies.
    Available via Nature Biotechnology supplementary materials (nbt.4061-S1.xlsx / Supplementary Table 1).
    """
    print("[*] Preparing dataset schemas for Kim et al. (2018) and DeWeirdt et al. (2021)...")
    # Schema definition for standardized guide evaluation
    schema = {
        "target_id": str,
        "full_34bp": str,      # 4bp upstream + 4bp PAM (TTTV) + 23bp protospacer + 3bp downstream
        "pam_seq": str,        # 4bp PAM
        "spacer_seq": str,     # 20-23bp guide sequence
        "efficiency_score": float, # Indel frequency or log2 fold change
        "cell_line": str,      # HEK293T, HCT116, etc.
        "dataset_source": str  # HT1, HT2, HT3, HEK-lenti, etc.
    }
    return schema

if __name__ == "__main__":
    setup_directories()
    fetch_kim_supplementary()
    print("[OK] Data acquisition module initialized.")
