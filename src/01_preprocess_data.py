"""
=============================================================================
CRISPR-Cas12a Explainable AI Project
Module: 01_preprocess_data.py
Author: Minh Tran (UWA / Independent Researcher)

Description:
  Reads the raw high-throughput screening data (Kim et al. 2018 Supplementary
  Table 1: HT1 dataset of 16,292 targets), performs quality control (QC),
  filters for valid PAM motifs (TTTV: TTTA, TTTC, TTTG), and outputs a
  standardized, clean dataset for downstream feature engineering.
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
RAW_DATA_PATH = os.path.join(BASE_DIR, "data", "raw", "Kim_2018_Supplementary_Table_1_HT1.xlsx")
PROCESSED_DIR = os.path.join(BASE_DIR, "data", "processed")
OUTPUT_CLEAN_CSV = os.path.join(PROCESSED_DIR, "clean_kim_2018_targets.csv")

def inspect_and_load_excel(file_path):
    print(f"[*] Loading raw dataset from:\n    {file_path}")
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"[!] Error: File not found at {file_path}")

    # Inspect first few rows to determine correct header row
    df_peek = pd.read_excel(file_path, sheet_name=0, nrows=2)
    has_title_row = (
        any("unnamed" in str(c).lower() for c in df_peek.columns) or 
        any("data set" in str(c).lower() for c in df_peek.columns)
    )
    
    header_idx = 1 if has_title_row else 0
    print(f"[+] Detected header row at index: {header_idx}")
    df = pd.read_excel(file_path, sheet_name=0, header=header_idx)
    print(f"[+] Loaded raw dataset with shape: {df.shape}")
    print(f"[+] Raw columns: {list(df.columns)}")
    return df

def clean_and_standardize_dataset(df):
    print("\n[*] Performing Data Quality Control (QC)...")
    
    # Identify sequence and indel columns flexibly
    target_seq_col = None
    indel_col = None
    
    for col in df.columns:
        c_lower = str(col).strip().lower()
        # Prefer explicit 34 bp column
        if ("34" in c_lower) and ("sequence" in c_lower or "target" in c_lower):
            target_seq_col = col
        elif ("target sequence" in c_lower or "sequence" in c_lower) and not target_seq_col:
            target_seq_col = col
            
        # Identify true indel frequency (%) column:
        # Strictly disqualify read counts, sequencing depth, total reads, and background controls
        is_count_col = any(term in c_lower for term in ["count", "read", "total", "depth"])
        is_background_col = "background" in c_lower
        is_freq_col = any(term in c_lower for term in ["freq", "%", "percent", "efficiency", "rate"])
        is_cpf1_col = any(term in c_lower for term in ["cpf1", "delivered", "ascas12a", "cas12a"])

        if is_freq_col and not is_count_col and not is_background_col:
            if is_cpf1_col:
                indel_col = col  # Primary choice: Cpf1-delivered indel frequency %
            elif not indel_col:
                indel_col = col

    # Fallback to inspecting content if column names are ambiguous
    if target_seq_col is None:
        for col in df.columns:
            sample = df[col].dropna().astype(str)
            if len(sample) > 0 and sample.str.len().mean() == 34.0:
                target_seq_col = col
                break

    if indel_col is None:
        for col in df.columns:
            c_lower = str(col).strip().lower()
            if not any(term in c_lower for term in ["count", "read", "total", "depth", "background"]):
                num = pd.to_numeric(df[col], errors="coerce").dropna()
                if len(num) > 0 and 20.0 < num.mean() < 70.0 and num.max() <= 100.0:
                    indel_col = col
                    break

    if target_seq_col is None:
        raise ValueError("Could not identify the 34 bp target sequence column in raw dataset.")
    if indel_col is None:
        raise ValueError("Could not identify the Cpf1-delivered indel frequency (%) column in raw dataset.")

    # Strict QC verification: indel frequency cannot exceed 100% (preventing read count contamination)
    num_check = pd.to_numeric(df[indel_col], errors="coerce").dropna()
    if num_check.max() > 100.0:
        raise ValueError(
            f"FATAL QC ERROR: Selected column '{indel_col}' contains values > 100 (max: {num_check.max()}). "
            "This is a read count column, not an indel percentage frequency column!"
        )

    print(f"[+] Identified sequence column: '{target_seq_col}'")
    print(f"[+] Identified indel frequency column: '{indel_col}' (mean: {num_check.mean():.2f}%, max: {num_check.max():.2f}%)")
    
    clean_df = pd.DataFrame()
    clean_df["raw_sequence"] = df[target_seq_col].astype(str).str.strip().str.upper()
    clean_df["indel_frequency"] = pd.to_numeric(df[indel_col], errors="coerce")
    
    # Filter valid 34 bp sequences
    initial_count = len(clean_df)
    clean_df = clean_df.dropna(subset=["raw_sequence", "indel_frequency"])
    
    # Ensure length is 34 bp (or extract 34bp if longer)
    clean_df = clean_df[clean_df["raw_sequence"].str.len() == 34]
    clean_df = clean_df[~clean_df["raw_sequence"].str.contains("[^ACGT]")]
    
    # Decompose sequence:
    # 4bp upstream [0:4], 4bp PAM [4:8], 23bp spacer [8:31], 3bp downstream [31:34]
    clean_df["upstream_4bp"] = clean_df["raw_sequence"].str.slice(0, 4)
    clean_df["pam_4bp"] = clean_df["raw_sequence"].str.slice(4, 8)
    clean_df["spacer_23bp"] = clean_df["raw_sequence"].str.slice(8, 31)
    clean_df["seed_8bp"] = clean_df["spacer_23bp"].str.slice(0, 8)
    clean_df["downstream_3bp"] = clean_df["raw_sequence"].str.slice(31, 34)
    
    # Flag valid canonical AsCas12a PAM: TTTV (where V is A, C, or G)
    clean_df["is_canonical_tttv"] = (
        clean_df["pam_4bp"].str.startswith("TTT") & 
        clean_df["pam_4bp"].str[3].isin(["A", "C", "G"])
    )
    
    # Binary classification label (consistent with literature: high activity >= 20% or top quartile)
    clean_df["is_high_efficiency"] = (clean_df["indel_frequency"] >= 20.0).astype(int)
    
    print(f"[OK] QC Complete:")
    print(f"    - Initial rows: {initial_count}")
    print(f"    - Valid 34bp rows retained: {len(clean_df)}")
    print(f"    - Canonical TTTV PAM count: {clean_df['is_canonical_tttv'].sum()} ({clean_df['is_canonical_tttv'].mean():.1%})")
    print(f"    - High-efficiency guides (indel >= 20%): {clean_df['is_high_efficiency'].sum()} ({clean_df['is_high_efficiency'].mean():.1%})")
    print(f"    - Indel frequency summary:\n{clean_df['indel_frequency'].describe()}")
    
    return clean_df

def main():
    os.makedirs(PROCESSED_DIR, exist_ok=True)
    try:
        raw_df = inspect_and_load_excel(RAW_DATA_PATH)
        clean_df = clean_and_standardize_dataset(raw_df)
        clean_df.to_csv(OUTPUT_CLEAN_CSV, index=False)
        print(f"\n[OK] Clean dataset successfully saved to:\n    {OUTPUT_CLEAN_CSV}")
    except Exception as e:
        print(f"[!] Processing failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
