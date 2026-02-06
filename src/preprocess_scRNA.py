import scanpy as sc
import sys
import importlib_metadata
import argparse
import warnings
import logging
import os

# Suppress annoying warnings while keeping functionality
warnings.filterwarnings("ignore", category=UserWarning, module="anndata")
logging.getLogger("scanpy").setLevel(logging.ERROR)

# Fix importlib issue in some environments
sys.modules["importlib.metadata"] = importlib_metadata


def read_samples(file_path):
    if not os.path.exists(file_path):
        raise FileNotFoundError(
            f"[ERROR] Samples file does not exist: {file_path}"
        )

    samples = {}

    with open(file_path, "r") as file:
        for line_num, line in enumerate(file, start=1):
            sample_id = line.strip()

            # Skip empty lines
            if not sample_id:
                continue

            # Reject CSV-like formatting
            if "," in sample_id:
                raise ValueError(
                    f"[ERROR] Invalid format in samples file at line {line_num}: "
                    f"'{sample_id}'. Expected ONE sample ID per line (no commas, no header)."
                )

            filename = f"{sample_id}_filtered_feature_bc_matrix.h5"
            samples[sample_id] = filename

    if len(samples) == 0:
        raise ValueError(
            "[ERROR] No valid samples found in samples file. "
            "File may be empty or incorrectly formatted."
        )

    return samples


def main():
    parser = argparse.ArgumentParser(
        description="Read multiple 10x samples, QC, and concatenate into one AnnData object"
    )

    parser.add_argument(
        "--project_name",
        required=True,
        help="Name of the output AnnData object (without .h5ad)"
    )

    parser.add_argument(
        "--samples",
        required=True,
        help="Text file with ONE sample ID per line"
    )

    args = parser.parse_args()

    project_name = args.project_name
    samples_file = args.samples

    print("[INFO] Starting preprocessing")
    print(f"[INFO] Project name: {project_name}")
    print(f"[INFO] Samples file: {samples_file}")

    # Read and validate samples
    samples = read_samples(samples_file)
    print(f"[INFO] Found {len(samples)} samples:")
    for s in samples:
        print(f"  - {s}")

    adatas = {}

    for sample_id, filename in samples.items():
        print(f"[INFO] Reading sample '{sample_id}' from {filename}")

        if not os.path.exists(filename):
            raise FileNotFoundError(
                f"[ERROR] Missing 10x file for sample '{sample_id}': {filename}"
            )

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata = sc.read_10x_h5(filename)

        # Fix duplicate names immediately
        adata.var_names_make_unique()
        adata.obs_names_make_unique()

        adata.obs["sample"] = sample_id
        adatas[sample_id] = adata

    if len(adatas) == 0:
        raise RuntimeError(
            "[ERROR] No AnnData objects were created. Cannot continue."
        )

    print("[INFO] Concatenating samples")
    combined_adata = sc.concat(
        adatas.values(),
        label="sample",
        keys=adatas.keys()
    )

    combined_adata.var_names_make_unique()
    combined_adata.obs_names_make_unique()

    # QC metrics
    combined_adata.var["mt"] = combined_adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(
        combined_adata,
        qc_vars=["mt"],
        inplace=True,
        log1p=True
    )

    # QC plots (before filtering)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pl.violin(
            combined_adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            save="_QC.png"
        )

    # Filtering
    combined_adata = combined_adata[
        (combined_adata.obs["n_genes_by_counts"] > 200) &
        (combined_adata.obs["n_genes_by_counts"] < 8000) &
        (combined_adata.obs["total_counts"] > 500) &
        (combined_adata.obs["total_counts"] < 30000) &
        (combined_adata.obs["pct_counts_mt"] < 20),
        :
    ].copy()

    sc.pp.filter_cells(combined_adata, min_genes=100)
    sc.pp.filter_genes(combined_adata, min_cells=3)

    # QC plots (after filtering)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc.pl.violin(
            combined_adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            save="_AfterQC.png"
        )

    # Save output
    output_file = f"{project_name}.h5ad"
    combined_adata.write(output_file)

    print(f"[SUCCESS] Created {output_file}")


if __name__ == "__main__":
    main()

