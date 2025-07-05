# 01_data_retrieval

## Overview
This module fetches all reported human acetylcholinesterase (AChE) inhibitors with IC₅₀ ≤ 100 nM from ChEMBL (target CHEMBL220), deduplicates to the best potency per compound, and exports a CSV for further analysis.

## Prerequisites
- Python 3.8 or later  
- Git (optional, for cloning the repository)

## Setup
```bash
# 1. Clone the repository (if not already done)
git clone https://github.com/yourusername/ache-inhibitor-pipeline.git
cd ache-inhibitor-pipeline/01_data_retrieval

# 2. (Optional) Create and activate a virtual environment
python -m venv venv
source venv/bin/activate    # Linux/macOS
venv\Scripts\activate       # Windows

# 3. Install dependencies
pip install -r requirements.txt

