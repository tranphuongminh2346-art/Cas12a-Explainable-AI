"""
Command-Line Interface (CLI) for cas12a-xai.

Enables on-target cleavage prediction, batch scoring, and explainable design critique
directly from the terminal.
"""

import sys
import os
import argparse
import json
import numpy as np
import pandas as pd

from . import __version__
from .predictor import Cas12aPredictor


def parse_args():
    parser = argparse.ArgumentParser(
        prog="cas12a-xai",
        description="cas12a-xai: Biophysically-Guided Explainable AI for CRISPR-Cas12a Cleavage Prediction."
    )
    parser.add_argument(
        "--seq", "-s",
        type=str,
        help="Input target DNA sequence (34 bp target site, 27 bp PAM+spacer, or 23 bp spacer)."
    )
    parser.add_argument(
        "--file", "-f",
        type=str,
        help="Path to input CSV file containing target sequences for batch prediction."
    )
    parser.add_argument(
        "--column", "-c",
        type=str,
        default="sequence",
        help="Name of the sequence column in input CSV (default: 'sequence')."
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        help="Path to save output predictions (CSV format)."
    )
    parser.add_argument(
        "--model", "-m",
        type=str,
        default=None,
        help="Optional path to a trained LightGBM model booster file."
    )
    parser.add_argument(
        "--explain", "-e",
        action="store_true",
        help="Print detailed biophysical explanation and design recommendations for the sequence."
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results in JSON format."
    )
    parser.add_argument(
        "--version", "-v",
        action="version",
        version=f"cas12a-xai v{__version__}"
    )
    return parser.parse_args()


def main():
    if hasattr(sys.stdout, 'reconfigure'):
        try:
            sys.stdout.reconfigure(encoding='utf-8')
        except Exception:
            pass

    args = parse_args()

    if not args.seq and not args.file:
        print("[-] Error: Please provide either a sequence via --seq or a file via --file.", file=sys.stderr)
        print("    Example (27 bp): cas12a-xai --seq TTTCGACCGTTAGGCTACGATCGGATC --explain", file=sys.stderr)
        print("    Example (34 bp): cas12a-xai --seq ACGGTTTCGACCGTTAGGCTACGATCGGATCGCC --explain", file=sys.stderr)
        sys.exit(1)

    predictor = Cas12aPredictor(model_path=args.model)

    # Single sequence mode
    if args.seq:
        seq = args.seq.strip()
        if args.explain:
            res = predictor.explain(seq)
            if args.json:
                print(json.dumps(res, indent=2))
            else:
                print("=" * 70)
                print(f"CRISPR-Cas12a Explainable AI Design Report")
                print("=" * 70)
                print(f"Target Sequence:       {res['target_sequence']}")
                print(f"Predicted Indel Eff:   {res['predicted_indel_percent']}% ({res['activity_tier']})")
                print(f"PAM Motif:             {res['parsed_domains']['pam_4bp']}")
                print(f"Spacer (23 bp):        {res['parsed_domains']['spacer_23bp']}")
                print("-" * 70)
                print("Biophysical Properties:")
                for k, v in res['biophysical_metrics'].items():
                    print(f"  - {k:<42}: {v}")
                
                if res['flags']:
                    print("-" * 70)
                    print("Attention Flags:")
                    for flag in res['flags']:
                        print(f"  [!] {flag}")

                print("-" * 70)
                print("Actionable Design Recommendations:")
                for rec in res['design_recommendations']:
                    print(f"  [+] {rec}")
                print("=" * 70)
        else:
            eff = float(predictor.predict([seq])[0])
            if args.json:
                print(json.dumps({"sequence": seq, "predicted_indel_percent": round(eff, 2)}))
            else:
                print(f"Predicted On-Target Indel Frequency: {eff:.2f}%")

    # Batch file mode
    elif args.file:
        if not os.path.exists(args.file):
            print(f"[-] Error: Input file '{args.file}' not found.", file=sys.stderr)
            sys.exit(1)

        df = pd.read_csv(args.file)
        if args.column not in df.columns:
            print(f"[-] Error: Column '{args.column}' not found in {args.file}.", file=sys.stderr)
            sys.exit(1)

        sequences = df[args.column].astype(str).tolist()
        print(f"[*] Processing {len(sequences):,} sequences...")
        preds = predictor.predict(sequences)
        df['predicted_efficiency'] = np.round(preds, 2)
        df['is_high_efficiency'] = (df['predicted_efficiency'] >= 20.0).astype(int)

        out_path = args.output if args.output else "cas12a_predictions.csv"
        df.to_csv(out_path, index=False)
        print(f"[OK] Saved predictions to '{out_path}' successfully!")


if __name__ == "__main__":
    main()
