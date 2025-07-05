#!/usr/bin/env python3
"""
Fetch human AChE inhibitors (CHEMBL220) with IC50 ≤100 nM from ChEMBL,
dedupe to the lowest IC50 per compound, and save as CSV.
"""

import logging
import pandas as pd
from chembl_webresource_client.new_client import new_client

# Configuration
TARGET_ID = "CHEMBL220"
IC50_THRESHOLD = 100  # nM
OUTPUT_CSV = "AChE_Inhibitors_IC50_leq100nM.csv"

# Logging
logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger()

def fetch_activities() -> pd.DataFrame:
    """Pull all IC50 ≤ threshold for the human AChE target."""
    client = new_client.activity
    activities = client.filter(
        target_chembl_id=TARGET_ID,
        standard_type="IC50",
        standard_units="nM",
        standard_value__lte=IC50_THRESHOLD,
        pchembl_value__isnull=False
    )
    
    records = []
    for act in activities:
        try:
            rel = act.get("standard_relation", "=")
            records.append({
                "CHEMBL_ID": act["molecule_chembl_id"],
                "SMILES": act["canonical_smiles"],
                "IC50_nM": float(act["standard_value"]),
                "Relation": rel,
                "Bioassay_ID": act["assay_chembl_id"],
                "Reference_ID": act.get("document_chembl_id", ""),
                "pChEMBL": float(act["pchembl_value"]) if act.get("pchembl_value") else None
            })
        except (KeyError, ValueError, TypeError):
            continue

    return pd.DataFrame(records)

def deduplicate_to_lowest_ic50(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only the lowest IC50 measurement per compound"""
    return (
        df.sort_values("IC50_nM")
        .groupby("CHEMBL_ID")
        .first()
        .reset_index()
    )

def main():
    log.info("Fetching AChE inhibitors with IC50 ≤100 nM...")
    df_raw = fetch_activities()
    
    if df_raw.empty:
        log.error("No records found matching criteria")
        return

    log.info("Retrieved %d activity records", len(df_raw))
    
    df_deduped = deduplicate_to_lowest_ic50(df_raw)
    log.info("Deduplicated to %d unique compounds", len(df_deduped))
    
    df_deduped.to_csv(OUTPUT_CSV, index=False)
    log.info("Saved results to %s", OUTPUT_CSV)

if __name__ == "__main__":
    main()
