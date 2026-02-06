import scanpy as sc
import pandas as pd
import sys
import importlib_metadata
import argparse
import os

# Fix for importlib metadata
sys.modules['importlib.metadata'] = importlib_metadata

# -----------------------------
# Parse arguments
# -----------------------------
parser = argparse.ArgumentParser(description="Export barcode annotations from AnnData")

parser.add_argument(
    '--project_name',
    required=True,
    help="Input AnnData file (.h5ad)"
)

parser.add_argument(
    '--barcodes',
    required=True,
    help="Output CSV file for barcode annotations"
)

args = parser.parse_args()

input_adata_file = args.project_name
output_barcode_csv = args.barcodes

# -----------------------------
# Load AnnData
# -----------------------------
adata = sc.read_h5ad(input_adata_file)

# -----------------------------
# Check obs columns
# -----------------------------
print(f"[INFO] Columns in adata.obs: {list(adata.obs.columns)}")

# -----------------------------
# Extract barcode annotations
# -----------------------------
# Adjust 'celltype' if your annotations are named differently
if 'celltype' not in adata.obs.columns:
    raise ValueError(
        "[ERROR] 'celltype' column not found in adata.obs. "
        "Available columns: " + ", ".join(adata.obs.columns)
    )

barcode_annotations = adata.obs[['celltype']].copy()
barcode_annotations.index.name = 'barcode'
barcode_annotations.reset_index(inplace=True)

# -----------------------------
# Save to CSV
# -----------------------------
barcode_annotations.to_csv(output_barcode_csv, index=False)
print(f"[SUCCESS] Saved barcode annotations to {output_barcode_csv}")

