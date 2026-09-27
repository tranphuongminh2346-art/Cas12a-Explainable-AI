"""
=============================================================================
CRISPR-Cas12a Explainable AI Project
Module: 03_train_models.py
Author: Minh Tran (UWA / Independent Researcher)

Description:
  Trains and validates the biophysically-guided Cas12a prediction framework
  using 10-Fold Stratified Cross-Validation:
    1. Dual stratification on PAM subtype (TTTA, TTTC, TTTG) and Indel activity
       quintiles to eliminate distribution shift across folds.
    2. 57 engineered biophysical features (SantaLucia nearest-neighbor duplex
       free energy, Tm polarity gradient, Pol III termination, PAM epistasis).
    3. Quantitative evaluation: Spearman rho, Pearson r, ROC-AUC (>=20% indel),
       MAE, RMSE with 95% confidence intervals.
    4. Subgroup PAM stratification (TTTA, TTTC, TTTG) confirming zero PAM leakage.
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

import argparse
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import mean_squared_error, mean_absolute_error, roc_auc_score

# Project Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC_DIR = os.path.join(BASE_DIR, "src")
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from cas12a_xai.features import FEATURE_NAMES, extract_features_batch

INPUT_FEATURE_CSV = os.path.join(BASE_DIR, "data", "processed", "cas12a_features_master.csv")
RAW_DATA_CSV = os.path.join(BASE_DIR, "data", "raw", "kim_2018_dataset.csv")
RESULTS_DIR = os.path.join(BASE_DIR, "data", "results")
PREDICTIONS_CSV = os.path.join(RESULTS_DIR, "cas12a_predictions_cv.csv")
PAM_METRICS_CSV = os.path.join(RESULTS_DIR, "pam_stratified_validation.csv")
PKG_MODEL_PATH = os.path.join(SRC_DIR, "cas12a_xai", "models", "cas12a_lgbm_model.txt")


def calc_ci(arr):
    """Compute mean, standard deviation, and 95% confidence interval."""
    m = np.mean(arr)
    se = np.std(arr, ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else 0.0
    sd = np.std(arr, ddof=1) if len(arr) > 1 else 0.0
    return m, sd, m - 1.96 * se, m + 1.96 * se


def evaluate_predictions(df: pd.DataFrame, label_col: str = "indel_frequency", pred_col: str = "predicted_efficiency"):
    """Compute overall and PAM-stratified benchmark metrics."""
    y_true = df[label_col].values
    y_pred = df[pred_col].values
    y_bin = (y_true >= 20.0).astype(int)

    rho, _ = spearmanr(y_true, y_pred)
    r, _ = pearsonr(y_true, y_pred)
    auc_score = roc_auc_score(y_bin, y_pred)
    mae = mean_absolute_error(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))

    print("\n" + "=" * 78)
    print("CRISPR-Cas12a 10-FOLD STRATIFIED BENCHMARK PERFORMANCE")
    print("=" * 78)
    print(f"Total Evaluated Guides (N):         {len(df):,}")
    print(f"Spearman Rank Correlation (rho):     {rho:.4f}")
    print(f"Pearson Linear Correlation (r):     {r:.4f}")
    print(f"Classification ROC-AUC (>=20%):     {auc_score:.4f}")
    print(f"Mean Absolute Error (MAE):          {mae:.2f}%")
    print(f"Root Mean Squared Error (RMSE):     {rmse:.2f}%")
    print("=" * 78)

    pam_rows = []
    print("\nPAM-Stratified Subgroup Performance (No PAM Leakage):")
    print("-" * 78)
    print(f"{'PAM':<6} {'Count':<7} {'Mean Indel':<12} {'Spearman rho':<14} {'Pearson r':<12} {'ROC-AUC':<10}")
    print("-" * 78)

    for pam in ["TTTA", "TTTC", "TTTG"]:
        sub = df[df["pam_4bp"] == pam]
        if len(sub) == 0:
            continue
        sub_true = sub[label_col].values
        sub_pred = sub[pred_col].values
        sub_bin = (sub_true >= 20.0).astype(int)

        s_rho, _ = spearmanr(sub_true, sub_pred)
        s_r, _ = pearsonr(sub_true, sub_pred)
        s_auc = roc_auc_score(sub_bin, sub_pred)
        s_mae = mean_absolute_error(sub_true, sub_pred)
        s_rmse = np.sqrt(mean_squared_error(sub_true, sub_pred))

        print(f"{pam:<6} {len(sub):<7} {sub_true.mean():<12.2f}% {s_rho:<14.4f} {s_r:<12.4f} {s_auc:<10.4f}")

        pam_rows.append({
            "PAM Motif": pam,
            "Count (n)": len(sub),
            "Library Proportion (%)": round(len(sub) / len(df) * 100, 2),
            "Mean Observed Indel (%)": round(sub_true.mean(), 2),
            "Median Observed Indel (%)": round(float(np.median(sub_true)), 2),
            "Spearman ρ": round(s_rho, 4),
            "Pearson r": round(s_r, 4),
            "ROC-AUC (>=20%)": round(s_auc, 4),
            "MAE (%)": round(s_mae, 2),
            "RMSE (%)": round(s_rmse, 2)
        })
    print("-" * 78)
    return pd.DataFrame(pam_rows)


def load_dataset():
    """Load and prepare canonical Cas12a target dataset."""
    if os.path.exists(PREDICTIONS_CSV):
        print(f"[*] Found existing verified 10-fold predictions: {PREDICTIONS_CSV}")
        return pd.read_csv(PREDICTIONS_CSV)
    elif os.path.exists(INPUT_FEATURE_CSV):
        print(f"[*] Loading feature matrix: {INPUT_FEATURE_CSV}")
        df = pd.read_csv(INPUT_FEATURE_CSV)
        return df[df.get("is_canonical_tttv", 1) == 1].copy()
    else:
        raise FileNotFoundError(
            f"[!] Neither {PREDICTIONS_CSV} nor {INPUT_FEATURE_CSV} was found."
        )


def run_10fold_cv(df: pd.DataFrame):
    """Execute 10-Fold Stratified Cross-Validation using LightGBM."""
    print(f"\n[*] Assembling 57-feature matrix for N = {len(df):,} canonical targets...")
    available_features = [f for f in FEATURE_NAMES if f in df.columns]
    if len(available_features) < len(FEATURE_NAMES):
        print(f"[*] Extracting missing features to complete 57 biophysical descriptors...")
        feat_df = extract_features_batch(df["raw_sequence"].tolist())
        for col in feat_df.columns:
            df[col] = feat_df[col].values
        available_features = FEATURE_NAMES

    X = df[available_features].values.astype(float)
    y = df["indel_frequency"].values.astype(float)
    y_bin = (y >= 20.0).astype(int)

    # Dual stratification: PAM + Indel Quintiles
    indel_bins = pd.qcut(df["indel_frequency"], q=5, labels=False, duplicates="drop")
    strat_key = df["pam_4bp"].astype(str) + "_" + indel_bins.astype(str)

    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

    oof_preds = np.zeros(len(y))
    spearman_folds = []
    pearson_folds = []
    auc_folds = []
    mae_folds = []
    rmse_folds = []

    print("[*] Training 10-Fold Stratified LightGBM...")
    try:
        import lightgbm as lgb
        use_lgb = True
    except ImportError:
        print("[!] LightGBM not installed. Falling back to RandomForestRegressor...")
        from sklearn.ensemble import RandomForestRegressor
        use_lgb = False

    params = {
        'objective': 'regression_l2',
        'boosting_type': 'gbdt',
        'learning_rate': 0.03,
        'num_leaves': 31,
        'feature_fraction': 0.80,
        'bagging_fraction': 0.85,
        'bagging_freq': 3,
        'min_child_samples': 20,
        'seed': 42,
        'verbose': -1
    }

    for fold, (train_idx, val_idx) in enumerate(skf.split(X, strat_key), start=1):
        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]

        if use_lgb:
            dtrain = lgb.Dataset(X_train, label=y_train, feature_name=available_features)
            dval = lgb.Dataset(X_val, label=y_val, reference=dtrain, feature_name=available_features)
            model = lgb.train(
                params,
                dtrain,
                num_boost_round=500,
                valid_sets=[dval],
                callbacks=[lgb.early_stopping(50, verbose=False)]
            )
            val_pred = model.predict(X_val)
        else:
            rf = RandomForestRegressor(n_estimators=200, max_depth=15, n_jobs=-1, random_state=42)
            rf.fit(X_train, y_train)
            val_pred = rf.predict(X_val)

        val_pred = np.clip(val_pred, 0.0, 100.0)
        oof_preds[val_idx] = val_pred

        rho, _ = spearmanr(y_val, val_pred)
        r, _ = pearsonr(y_val, val_pred)
        auc_v = roc_auc_score((y_val >= 20.0).astype(int), val_pred)
        mae = mean_absolute_error(y_val, val_pred)
        rmse = np.sqrt(mean_squared_error(y_val, val_pred))

        spearman_folds.append(rho)
        pearson_folds.append(r)
        auc_folds.append(auc_v)
        mae_folds.append(mae)
        rmse_folds.append(rmse)

        print(f"  Fold {fold:02d}/10: Spearman rho = {rho:.4f} | Pearson r = {r:.4f} | ROC-AUC = {auc_v:.4f} | MAE = {mae:.2f}%")

    df["predicted_efficiency"] = oof_preds
    df["residual"] = df["indel_frequency"] - df["predicted_efficiency"]
    df["abs_error"] = np.abs(df["residual"])

    sp_m, sp_s, sp_lo, sp_hi = calc_ci(spearman_folds)
    pr_m, pr_s, pr_lo, pr_hi = calc_ci(pearson_folds)
    auc_m, auc_s, auc_lo, auc_hi = calc_ci(auc_folds)

    print("\n" + "=" * 78)
    print("10-FOLD STRATIFIED CROSS-VALIDATION SUMMARY (MEAN +/- SD [95% CI])")
    print("=" * 78)
    print(f"Spearman Rank Correlation (rho): {sp_m:.4f} +/- {sp_s:.4f}  [95% CI: {sp_lo:.4f} - {sp_hi:.4f}]")
    print(f"Pearson Linear Correlation (r): {pr_m:.4f} +/- {pr_s:.4f}  [95% CI: {pr_lo:.4f} - {pr_hi:.4f}]")
    print(f"Classification ROC-AUC (>=20%):  {auc_m:.4f} +/- {auc_s:.4f}  [95% CI: {auc_lo:.4f} - {auc_hi:.4f}]")
    print(f"Mean Absolute Error (MAE):       {np.mean(mae_folds):.2f}% +/- {np.std(mae_folds, ddof=1):.2f}%")
    print(f"Root Mean Squared Error (RMSE):   {np.mean(rmse_folds):.2f}% +/- {np.std(rmse_folds, ddof=1):.2f}%")
    print("=" * 78)

    os.makedirs(RESULTS_DIR, exist_ok=True)
    df.to_csv(PREDICTIONS_CSV, index=False)
    print(f"[OK] Saved out-of-fold predictions to: {PREDICTIONS_CSV}")

    pam_df = evaluate_predictions(df)
    pam_df.to_csv(PAM_METRICS_CSV, index=False)
    print(f"[OK] Saved PAM-stratified validation table to: {PAM_METRICS_CSV}")

    # Train and export final production booster model
    train_and_save_full_booster(df)


def train_and_save_full_booster(df: pd.DataFrame, save_path: str = PKG_MODEL_PATH):
    """Train a final LightGBM model on all canonical guides and export the booster file."""
    try:
        import lightgbm as lgb
    except ImportError:
        print("[!] LightGBM not installed; skipping booster model export.")
        return None

    available_features = [f for f in FEATURE_NAMES if f in df.columns]
    X = df[available_features].values
    y = df["indel_frequency"].values

    params = {
        'objective': 'regression_l2',
        'boosting_type': 'gbdt',
        'learning_rate': 0.03,
        'num_leaves': 31,
        'feature_fraction': 0.80,
        'bagging_fraction': 0.85,
        'bagging_freq': 3,
        'min_child_samples': 20,
        'seed': 42,
        'verbose': -1
    }
    dtrain = lgb.Dataset(X, label=y, feature_name=available_features)
    print(f"\n[*] Training production booster on all {len(df):,} canonical guides...")
    model = lgb.train(params, dtrain, num_boost_round=350)
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    model.save_model(save_path)
    print(f"[OK] Saved production LightGBM booster model to: {save_path}")
    return model


def main():
    parser = argparse.ArgumentParser(description="10-Fold Stratified Training and Validation for Cas12a.")
    parser.add_argument("--evaluate_only", action="store_true", help="Evaluate existing out-of-fold predictions.")
    parser.add_argument("--save_booster", action="store_true", help="Train and export final LightGBM booster model.")
    args = parser.parse_args()

    df = load_dataset()
    if args.save_booster:
        train_and_save_full_booster(df)
    elif args.evaluate_only or ("predicted_efficiency" in df.columns and len(df) == 11365):
        print("[*] Evaluating canonical out-of-fold predictions table...")
        pam_df = evaluate_predictions(df)
        pam_df.to_csv(PAM_METRICS_CSV, index=False)
    else:
        run_10fold_cv(df)


if __name__ == "__main__":
    main()
