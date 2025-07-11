#!/usr/bin/env python3
"""
admet_filter.py

Apply strict ADMET + Lipinski Rule of 5 filters (with a TPSA ≥ 20 Å² floor)
AND loose pharmacokinetic property filters to the deduplicated ADMETlab 3.0 dataset.
"""

import pandas as pd
import os
import sys

INPUT_CSV  = "admet_full_dedup.csv"
OUTPUT_CSV = "admet_filtered.csv"

def main():
    if not os.path.isfile(INPUT_CSV):
        print(f"Error: `{INPUT_CSV}` not found.")
        sys.exit(1)

    df = pd.read_csv(INPUT_CSV)

    strict_filters = (
        (df['logS']   >= -4.5) &
        (df['logP']   <=  5.0) &
        (df['TPSA']   >= 20.0) &   # require at least 20 Å² polar surface
        (df['TPSA']   <=100.0) &
        (df['nRot']   <= 10  ) &
        (df['BBB']    >=  0.70) &
        (df['MW']     <=500.0) &
        (df['nHD']    <=  5  ) &
        (df['nHA']    <= 10   )
    )
    df_strict = df[strict_filters]

    loose_filters = (
        (df_strict['logVDss']   <=  5) &
        (df_strict['cl-plasma'] <= 50) &
        (df_strict['t0.5']      >=  0.3) &
        (df_strict['PPB']       <= 99.5) &
        (df_strict['Fsp3']      >=  0.1)
    )
    df_final = df_strict[loose_filters]

    df_final.to_csv(OUTPUT_CSV, index=False)

    print(f"Total molecules:           {len(df)}")
    print(f"After strict filtering:    {len(df_strict)}")
    print(f"After loose filtering:     {len(df_final)}")
    print(f"Results saved to `{OUTPUT_CSV}`")

if __name__ == "__main__":
    main()
