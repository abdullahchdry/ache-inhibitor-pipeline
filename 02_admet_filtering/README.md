# 02_admet_filtering

## Overview
This module applies strict ADMET + Lipinski filtering (no toxicity) to the deduplicated ADMETlab 3.0 export. It reads `admet_full_dedup.csv`, filters on:
- log S ≥ –4.5  
- log P ≤ 5.0  
- TPSA ≤ 100 Å²  
- nRot ≤ 10  
- BBB probability ≥ 0.70  
- MW ≤ 500 Da  
- H-bond donors ≤ 5  
- H-bond acceptors ≤ 10  

and writes the result to `admet_filtered.csv` for downstream UMAP/clustering.

## Prerequisites
- Python 3.8 or later  
- pandas  

## Setup & Run
```bash
# 1. Navigate to this module’s folder
cd ache-inhibitor-pipeline/02_admet_filtering

# 2. (Optional) Activate your virtual environment
source ../venv/bin/activate    # Linux/macOS
../venv\Scripts\activate       # Windows

# 3. Install dependencies (if you haven’t already)
pip install -r ../requirements.txt

# 4. Ensure your input CSV is present:
#    ../01_data_retrieval/admet_full_dedup.csv → copy or symlink here as admet_full_dedup.csv

# 5. Run the filter script
python admet_filter.py

# 6. Check output
#    admet_filtered.csv will be created in this directory
