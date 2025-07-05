#!/usr/bin/env python3
"""
Fetch human AChE inhibitors (CHEMBL220) with IC50 ≤100 nM from ChEMBL,
dedupe to the lowest IC50 per compound, and save as CSV.
"""

import logging
import pandas as pd
from chembl_webresource_client.new_client import new_client
import os

# Configuration
TARGET_ID       = "CHEMBL220"
IC50_THRESHOLD  = 100  # nM
OUTPUT_CSV      = "AChE_Inhibitors_IC50_leq100nM.csv"
FIELDS          = [
    "molecule_chembl_id", "canonical_smiles", "standard_value",
    "standard_relation", "assay_chembl_id", "document_chembl_id",
    "pchembl_value"
]

# Logging
logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger()

def fetch_activities() -> pd.DataFrame:
    """Pull all IC50 ≤ threshold for the human AChE target."""
    client = new_client.activity
    iterator = client.filter(
        target_chembl_id=TARGET_ID,
        standard_type="IC50",
        standard_units="nM",
        standard_value__lte=IC50_THRESHOLD,
        pchembl_value__isnull=False
    ).iterfilter()

    records = []
    for act in iterator:
        rel = act.get("standard_relation", "=")
        try:
            records.append({
                "CHEMBL_ID":      act["molecule_chembl_id"],
                "SMILES":         act.get("canonical_smiles"),
                "IC50_nM":        float(act["standard_value"]),
                "Relation":       rel,
                "Bioassay_ID":    act.get("assay_chembl_id"),
                "Reference_ID":   act.get("document_chembl_id", ""),
                "pChEMBL":        float(act["pchembl_value"])
            })
        except (KeyError, ValueError, TypeError):
            continue

    return pd.DataFrame(records)

def dedupe(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only '=' or '<=' relations, then the lowest IC50 per compound."""
    df = df[df["Relation"].isin(["=", "<="])]
    df = df.sort_values("IC50_nM").drop_duplicates("CHEMBL_ID", keep="first")
    return df

def main():
    log.info("Fetching raw activities...")
    df_raw = fetch_activities()
    if df_raw.empty:
        log.error("No records found. Check your filters.")
        return

    log.info("Deduplicating to best IC50 per compound...")
    df_best = dedupe(df_raw)

    df_best.to_csv(OUTPUT_CSV, index=False)
    log.info("Saved %d unique compounds to %s",
             len(df_best), OUTPUT_CSV)

if __name__ == "__main__":
    main()
