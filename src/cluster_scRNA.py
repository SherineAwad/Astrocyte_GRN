import scanpy as sc
import sys
import importlib_metadata
import matplotlib.pyplot as plt
import argparse
import os

# Fix for importlib metadata
sys.modules['importlib.metadata'] = importlib_metadata

# -----------------------------
# Parse arguments
# -----------------------------
parser = argparse.ArgumentParser(description="Cluster and plot scRNA-seq AnnData object")

parser.add_argument(
    '--project_name',
    required=True,
    help="Input AnnData file (.h5ad)"
)

parser.add_argument(
    '--markers',
    required=True,
    help="Text file with one marker gene per line"
)

parser.add_argument(
    '--output',
    required=True,
    help="Output AnnData filename (e.g. clustered_astrocyte.h5ad)"
)

parser.add_argument(
    '--prefix',
    required=True,
    help="Prefix for all figure files (e.g. astrocyte)"
)

args = parser.parse_args()

input_adata_file = args.project_name
markers_file = args.markers
output_adata_file = args.output
fig_prefix = args.prefix

print(f"[INFO] Input AnnData: {input_adata_file}")
print(f"[INFO] Markers file: {markers_file}")
print(f"[INFO] Output AnnData: {output_adata_file}")
print(f"[INFO] Figure prefix: {fig_prefix}")

# Base name for figure outputs
base_name = os.path.splitext(os.path.basename(fig_prefix))[0]

# Create figures directory if it doesn't exist
os.makedirs("figures", exist_ok=True)

# -----------------------------
# Load AnnData
# -----------------------------
combined_adata = sc.read(input_adata_file)

# -----------------------------
# CRITICAL FIX: Preserve raw counts for SCENIC+
# -----------------------------
# Store the raw counts before normalization
combined_adata.raw = combined_adata
print("✅ Raw counts preserved in adata.raw for SCENIC+ compatibility")

# -----------------------------
# Rename samples before anything else
# -----------------------------
rename_dict = {"TH1": "Control", "TH2": "OE"}
if "sample" in combined_adata.obs.columns:
    combined_adata.obs["sample"] = combined_adata.obs["sample"].replace(rename_dict)
    print("✅ Renamed samples:", combined_adata.obs["sample"].unique())
else:
    print("⚠️ Column 'sample' not found. Available columns:", combined_adata.obs.columns.tolist())

# -----------------------------
# Preprocessing & clustering
# -----------------------------
sc.pp.normalize_total(combined_adata, target_sum=1e4)
sc.pp.log1p(combined_adata)
combined_adata.raw = combined_adata
sc.pp.scale(combined_adata)

sc.tl.pca(combined_adata, svd_solver="arpack")
sc.pp.neighbors(combined_adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(combined_adata)

sc.tl.leiden(combined_adata, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)

# Save cluster plot
cluster_fig_name = f"_{fig_prefix}_Clusters.png"
sc.pl.umap(combined_adata, color=["leiden"], legend_loc="on data", save=cluster_fig_name)
sample_fig_name = f"_{fig_prefix}_Samples.png"
sc.pl.umap(combined_adata, color="sample", legend_loc="on data", save=sample_fig_name)

# -----------------------------
# Marker gene plotting
# -----------------------------
with open(markers_file) as f:
    marker_genes = [line.strip() for line in f]

figurename = f"figures/{fig_prefix}_Dotplot.png"
combined_adata.obs["leiden"] = combined_adata.obs["leiden"].astype("category")
marker_genes_present = [g for g in marker_genes if g in combined_adata.var_names]

print(f"[INFO] Markers found in dataset: {marker_genes_present}")

fig = sc.pl.dotplot(
    combined_adata,
    var_names=marker_genes_present,
    groupby="leiden",
    standard_scale="var",
    show=False,
    return_fig=True
)
# Save figure
fig.savefig(figurename, dpi=600, bbox_inches="tight")

# -----------------------------
# Individual gene UMAP plots
# -----------------------------
for gene in marker_genes:
    if gene in combined_adata.var_names:
        sc.pl.umap(
            combined_adata,
            color=gene,
            title=gene,
            save=f"_{fig_prefix}_{gene}.png"
        )
    else:
        print(f"⚠️ Skipping {gene}: not in adata.var_names")

# -----------------------------
# Save clustered object
# -----------------------------
combined_adata.obs_names_make_unique()
combined_adata.write(output_adata_file, compression="gzip")
print(f"✅ Saved clustered AnnData to {output_adata_file}")
print(f"✅ Raw counts preserved for SCENIC+ compatibility")
