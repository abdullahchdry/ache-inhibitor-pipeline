#!/usr/bin/env python3
"""
admet_filter.py

Apply strict ADMET + Lipinski Rule of 5 filters (no toxicity) 
AND loose pharmacokinetic property filters to the deduplicated ADMETlab 3.0 dataset.

Requirements:
- Place your ADMETlab 3.0 export file named `admet_full_dedup.csv`
  in the same directory as this script. It must include columns:
  'logS', 'logP', 'TPSA', 'nRot', 'BBB', 'MW', 'nHD', 'nHA',
  'logVDss', 'cl-plasma', 't0.5', 'PPB', 'Fsp3'

Usage:
    python admet_filter.py

Output:
    Writes filtered results to `admet_filtered.csv`.
"""

import pandas as pd
import os
import sys

# File names
INPUT_CSV  = "admet_full_dedup.csv"
OUTPUT_CSV = "admet_filtered.csv"

def main():
    # Verify input file exists
    if not os.path.isfile(INPUT_CSV):
        print(f"Error: `{INPUT_CSV}` not found. Please place your ADMETlab 3.0 export in this directory with that name.")
        sys.exit(1)

    # Load the ADMETlab export
    df = pd.read_csv(INPUT_CSV)

    # Apply strict ADMET + Lipinski filtering
    strict_filters = (
        (df['logS']  >= -4.5) &       # Solubility
        (df['logP']  <=  5.0) &       # Lipophilicity (Lipinski)
        (df['TPSA']  <= 100.0) &      # Polarity
        (df['nRot']  <= 10) &         # Flexibility
        (df['BBB']   >= 0.70) &       # CNS penetration
        (df['MW']    <= 500.0) &      # Molecular weight
        (df['nHD']   <= 5) &          # H-bond donors
        (df['nHA']   <= 10)           # H-bond acceptors
    )

    df_strict = df[strict_filters]

    # Apply loose pharmacokinetic filters (soft limits)
    loose_filters = (
        (df_strict['logVDss'] <= 5) &          # Distribution volume
        (df_strict['cl-plasma'] <= 50) &       # Clearance
        (df_strict['t0.5'] >= 0.3) &           # Half-life
        (df_strict['PPB'] <= 99.5) &           # Plasma protein binding
        (df_strict['Fsp3'] >= 0.1)             # 3D complexity
    )

    df_final = df_strict[loose_filters]

    # Save filtered results
    df_final.to_csv(OUTPUT_CSV, index=False)

    print(f"Total molecules in dataset: {len(df)}")
    print(f"After strict filtering: {len(df_strict)} compounds")
    print(f"After loose additional filtering: {len(df_final)} compounds")
    print(f"Results saved to `{OUTPUT_CSV}`")

if __name__ == "__main__":
    main()
