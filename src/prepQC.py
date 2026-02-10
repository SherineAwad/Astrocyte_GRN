#!/usr/bin/env python3
import os
import argparse
import pickle

def read_samples_map(samples_map_path):
    fragments_dict = {}
    with open(samples_map_path) as f:
        for line in f:
            line = line.strip()
            if not line or ":" not in line:
                continue
            prefix, sample_name = [x.strip() for x in line.split(":", 1)]
            fragment_file = f"{prefix}_atac_fragments.tsv.gz"
            fragments_dict[sample_name] = fragment_file
    return fragments_dict

def main():
    parser = argparse.ArgumentParser(description="Generate pycisTopic QC command file for multiple samples")
    parser.add_argument("--out_dir", required=True, help="Base output directory")
    parser.add_argument("--consensus_dir", required=True, help="Directory containing consensus peak calling bed files")
    parser.add_argument("--tss_bed", required=True, help="Path to TSS BED file (e.g., tss_mm10.bed)")
    parser.add_argument("--samples_map", required=True, help="Path to samples_map.txt file")
    parser.add_argument("--qc_commands_filename", default="pycistopic_qc_commands.txt",
                        help="Filename to save QC commands")
    parser.add_argument("--min_fragments_per_cb", type=int, default=10,
                        help="Minimum number of fragments per cell barcode to keep (to avoid KDE errors). Default: 10")
    args = parser.parse_args()

    fragments_dict = read_samples_map(args.samples_map)

    fragments_dict_path = os.path.join(args.out_dir, "fragments_dict.pkl")
    os.makedirs(args.out_dir, exist_ok=True)
    with open(fragments_dict_path, "wb") as f:
        pickle.dump(fragments_dict, f)
    print(f"Saved fragments_dict to {fragments_dict_path}")

    consensus_regions_dir = args.consensus_dir
    if not os.path.exists(consensus_regions_dir):
        raise FileNotFoundError(f"Consensus regions directory not found: {consensus_regions_dir}")

    bed_files = [f for f in os.listdir(consensus_regions_dir) if f.endswith(".bed")]
    if len(bed_files) == 0:
        raise FileNotFoundError(f"No BED files found in {consensus_regions_dir}")

    bed_dict = {}
    for bed_file in bed_files:
        sample_name_from_bed = os.path.splitext(bed_file)[0]
        bed_dict[sample_name_from_bed] = os.path.join(consensus_regions_dir, bed_file)

    tss_bed_filename = args.tss_bed
    if not os.path.exists(tss_bed_filename):
        raise FileNotFoundError(f"TSS BED file not found: {tss_bed_filename}")

    qc_output_dir = os.path.join(args.out_dir, "QC")
    os.makedirs(qc_output_dir, exist_ok=True)

    with open(args.qc_commands_filename, "w") as fh:
        for sample, fragment_filename in fragments_dict.items():
            sample_qc_dir = os.path.join(qc_output_dir, sample)
            os.makedirs(sample_qc_dir, exist_ok=True)

            if sample not in bed_dict:
                raise FileNotFoundError(f"No consensus BED file matching sample '{sample}' found in {consensus_regions_dir}")

            regions_bed_filename = bed_dict[sample]

            # Compose full command as a single line string with all options together
            command = (
                "pycistopic qc "
                f"--fragments {fragment_filename} "
                f"--regions {regions_bed_filename} "
                f"--tss {tss_bed_filename} "
                f"--output {sample_qc_dir} "
                f"--min_fragments_per_cb {args.min_fragments_per_cb}"
            )
            fh.write(command + "\n")  # write as one line per sample command

    print(f"QC commands written to {args.qc_commands_filename}")
    print("You can now run them in the terminal, e.g.:")
    print(f"  bash {args.qc_commands_filename}")

if __name__ == "__main__":
    main()

