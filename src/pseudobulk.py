import os
import argparse
import pandas as pd
import pycisTopic
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk

def read_samples_map(file_path):
    """
    Reads sample map file with lines like:
      13005_TH2: KO1
      13784-TH1: Control
      13784-TH2: KO2
    Returns dict {sample_name: fragments_file}
    Assumes fragments files are named "<prefix>_atac_fragments.tsv.gz"
    """
    fragments_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if ':' not in line:
                print(f"Skipping malformed line: {line}")
                continue
            prefix, sample_name = line.split(':', 1)
            prefix = prefix.strip()
            sample_name = sample_name.strip()
            fragments_file = f"{prefix}_atac_fragments.tsv.gz"
            fragments_dict[sample_name] = fragments_file
    return fragments_dict

def main():
    parser = argparse.ArgumentParser(description="scenic+ pseudobulk peak calling with flexible samples")
    parser.add_argument('--samples_map', required=True, help="Path to sample map file")
    parser.add_argument('--out_dir', default="scenicOuts", help="Output directory")
    parser.add_argument('--chromsizes', default="https://hgdownload.cse.ucsc.edu/goldenpath/mm10/bigZips/mm10.chrom.sizes", help="Chromsizes file or URL")
    parser.add_argument('--n_cpu', type=int, default=4, help="Number of CPUs to use")
    parser.add_argument('--temp_dir', default="/tmp", help="Temporary directory for processing")
    args = parser.parse_args()

    print("pycisTopic version:", pycisTopic.__version__)

    # Step 1: Read fragments dictionary from sample map file
    fragments_dict = read_samples_map(args.samples_map)
    print("Fragments dictionary loaded:")
    for sample, frag_file in fragments_dict.items():
        print(f"  Sample: {sample} --> Fragment file: {frag_file}")

    # Step 2: Load chromsizes
    chromsizes = pd.read_table(
        args.chromsizes,
        header=None,
        names=["Chromosome", "End"]
    )
    chromsizes.insert(1, "Start", 0)
    print("\nChromsizes head:")
    print(chromsizes.head())

    # Step 3: Create output directories
    consensus_dir = os.path.join(args.out_dir, "consensus_peak_calling")
    bed_dir = os.path.join(consensus_dir, "pseudobulk_bed_files")
    bw_dir = os.path.join(consensus_dir, "pseudobulk_bw_files")

    os.makedirs(consensus_dir, exist_ok=True)
    os.makedirs(bed_dir, exist_ok=True)
    os.makedirs(bw_dir, exist_ok=True)

    # Step 4: Create cell data from fragments files
    cell_data_list = []

    for sample, fragments_file in fragments_dict.items():
        print(f"\nReading barcodes from {fragments_file} for sample {sample}...")

        if not os.path.exists(fragments_file):
            print(f"  ERROR: File not found: {fragments_file}")
            continue

        try:
            # Read fragments file - SKIP HEADER LINES starting with #
            fragments_df = pd.read_csv(
                fragments_file,
                sep='\t',
                header=None,
                comment='#',
                names=['chrom', 'start', 'end', 'barcode', 'count'],
                usecols=[0, 1, 2, 3, 4]
            )

            print(f"  Data shape (after skipping headers): {fragments_df.shape}")
            print(f"  First few barcodes: {fragments_df['barcode'].head().tolist()}")

            # Get unique barcodes
            barcodes = fragments_df['barcode'].unique()
            print(f"  Found {len(barcodes)} unique barcodes")

            # Create DataFrame for this sample
            sample_df = pd.DataFrame({
                'barcode': barcodes.astype(str),
                'sample': str(sample),
                'celltype': str(sample)
            })

            cell_data_list.append(sample_df)
            print(f"  Successfully processed {len(barcodes)} barcodes for sample {sample}")

        except Exception as e:
            print(f"  Error processing {fragments_file}: {e}")
            # Try alternative approach if the first fails
            try:
                print("  Trying alternative reading method...")
                fragments_df = pd.read_csv(
                    fragments_file,
                    sep='\t',
                    header=None,
                    comment='#',
                    usecols=[3]  # Only barcode column
                )
                fragments_df.columns = ['barcode']
                barcodes = fragments_df['barcode'].unique()
                print(f"  Alternative method found {len(barcodes)} barcodes")

                sample_df = pd.DataFrame({
                    'barcode': barcodes.astype(str),
                    'sample': str(sample),
                    'celltype': str(sample)
                })
                cell_data_list.append(sample_df)

            except Exception as e2:
                print(f"  Alternative method also failed: {e2}")

    if len(cell_data_list) == 0:
        print("\nERROR: No barcode data found in any fragment file!")
        exit(1)

    # Combine all samples
    cell_data = pd.concat(cell_data_list, ignore_index=True)

    # Create simplified version with only required columns
    cell_data_simple = cell_data[['barcode', 'sample', 'celltype']].copy()

    # Ensure all columns are strings
    for col in ['barcode', 'sample', 'celltype']:
        cell_data_simple[col] = cell_data_simple[col].astype(str)

    print(f"\nFinal cell data summary:")
    print(f"Total barcodes across all samples: {len(cell_data_simple)}")
    print("Sample distribution:")
    print(cell_data_simple['sample'].value_counts())
    print("\nFirst few rows:")
    print(cell_data_simple.head())

    # Step 5: Run pseudobulk export
    print("\nStarting pseudobulk export...")

    try:
        bw_paths, bed_paths = export_pseudobulk(
            input_data=cell_data_simple,
            variable="celltype",
            sample_id_col="sample",
            chromsizes=chromsizes,
            bed_path=bed_dir,
            bigwig_path=bw_dir,
            path_to_fragments=fragments_dict,
            n_cpu=args.n_cpu,
            normalize_bigwig=True,
            temp_dir=args.temp_dir
        )

        print("Pseudobulk export completed successfully!")

        # Save paths
        with open(os.path.join(consensus_dir, "bw_paths.tsv"), "w") as f:
            for v in bw_paths:
                f.write(f"{v}\t{bw_paths[v]}\n")

        with open(os.path.join(consensus_dir, "bed_paths.tsv"), "w") as f:
            for v in bed_paths:
                f.write(f"{v}\t{bed_paths[v]}\n")

        print(f"\nOutput directories:")
        print(f"  - BED files: {bed_dir}")
        print(f"  - BigWig files: {bw_dir}")
        print(f"  - Path lists: {consensus_dir}")

        print("\nGenerated files summary:")
        for celltype_sample, bed_path in bed_paths.items():
            if os.path.exists(bed_path):
                file_size = os.path.getsize(bed_path) / (1024 * 1024)
                print(f"  - {celltype_sample}: {bed_path} ({file_size:.2f} MB)")

    except Exception as e:
        print(f"Error during pseudobulk export: {e}")

        # Final attempt with minimal parameters
        print("\nTrying with minimal parameters...")
        try:
            bw_paths, bed_paths = export_pseudobulk(
                input_data=cell_data_simple,
                variable="celltype",
                sample_id_col="sample",
                chromsizes=chromsizes,
                bed_path=bed_dir,
                bigwig_path=bw_dir,
                path_to_fragments=fragments_dict,
                n_cpu=2
            )
            print("Success with minimal parameters!")

            # Save paths
            with open(os.path.join(consensus_dir, "bw_paths.tsv"), "w") as f:
                for v in bw_paths:
                    f.write(f"{v}\t{bw_paths[v]}\n")

            with open(os.path.join(consensus_dir, "bed_paths.tsv"), "w") as f:
                for v in bed_paths:
                    f.write(f"{v}\t{bed_paths[v]}\n")

        except Exception as e2:
            print(f"Final error: {e2}")

if __name__ == "__main__":
    main()

