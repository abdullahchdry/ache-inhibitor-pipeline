#!/usr/bin/env python3
"""
Fetch human AChE inhibitors (CHEMBL220) with IC50 ≤100 nM from ChEMBL,
override CHEMBL2448138 SMILES with the given parent string,
strip salts to the parent form, dedupe to the lowest IC50 per compound,
and save as CSV.

Note: CHEMBL2448138 lacked a valid SMILES in ChEMBL, so we hard-coded
its canonical parent SMILES at the end of this script.
"""

import logging
import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem

# Configuration
TARGET_ID   = "CHEMBL220"
IC50_THRESH = 100   # nM
OUTPUT_CSV  = "AChE_Inhibitors_IC50_leq100nM_parent_override.csv"

# Hard-coded override for CHEMBL2448138
OVERRIDE_ID     = "CHEMBL2448138"
OVERRIDE_SMILES = "COc1cc2c(cc1OC)C(=O)C(CC1CCN(CCNc3c4c(nc5ccccc35)CCCC4)CC1)C2"

# Logging
logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger()

def strip_salt(smiles: str) -> str:
    """Keep only the largest fragment by heavy-atom count, return canonical SMILES."""
    if not smiles:
        return ""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    parent = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    return Chem.MolToSmiles(parent, canonical=True)

def fetch_activities() -> pd.DataFrame:
    """Pull all IC50 ≤ threshold for the human AChE target, overriding as needed."""
    client     = new_client.activity
    activities = client.filter(
        target_chembl_id=TARGET_ID,
        standard_type="IC50",
        standard_units="nM",
        standard_value__lte=IC50_THRESH,
        pchembl_value__isnull=False
    )

    records = []
    for act in activities:
        chembl_id = act.get("molecule_chembl_id")
        raw_smiles = act.get("canonical_smiles") or ""
        # override CHEMBL2448138
        if chembl_id == OVERRIDE_ID:
            smiles = OVERRIDE_SMILES
        else:
            smiles = strip_salt(raw_smiles)

        try:
            ic50 = float(act["standard_value"])
        except (TypeError, ValueError):
            continue

        records.append({
            "CHEMBL_ID":    chembl_id,
            "SMILES":       smiles,
            "IC50_nM":      ic50,
            "Relation":     act.get("standard_relation", ""),
            "Bioassay_ID":  act.get("assay_chembl_id",""),
            "Reference_ID": act.get("document_chembl_id",""),
            "pChEMBL":      float(act["pchembl_value"]) if act.get("pchembl_value") else None
        })

    return pd.DataFrame(records)

def dedupe_ic50(df: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate to lowest IC50 per CHEMBL_ID."""
    return (
        df.sort_values("IC50_nM", ascending=True)
          .groupby("CHEMBL_ID", as_index=False)
          .first()
    )

def dedupe_smiles(df: pd.DataFrame) -> pd.DataFrame:
    """Drop any remaining exact SMILES duplicates, keeping the first."""
    return df.drop_duplicates(subset="SMILES", keep="first").reset_index(drop=True)

def main():
    log.info("Fetching AChE inhibitors with IC50 ≤ %d nM...", IC50_THRESH)
    df_raw = fetch_activities()
    log.info("Retrieved %d activity records", len(df_raw))

    df_by_ic50 = dedupe_ic50(df_raw)
    log.info("After IC50 dedupe: %d unique CHEMBL IDs", len(df_by_ic50))

    df_final = dedupe_smiles(df_by_ic50)
    log.info("After SMILES dedupe: %d unique parent molecules", len(df_final))

    df_final.to_csv(OUTPUT_CSV, index=False)
    log.info("Saved cleaned results to %s", OUTPUT_CSV)

    # Final note about the hard-coded override
    log.info("Note: CHEMBL2448138 parent SMILES was hard-coded as it was missing in ChEMBL.")

if __name__ == "__main__":
    main()
