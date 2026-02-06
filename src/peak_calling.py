#!/usr/bin/env python3
import os
import subprocess
import argparse

def read_samples_map(file_path):
    """
    Reads sample map file with lines like:
    13005_TH2: KO1
    13784-TH1: Control
    13784-TH2: KO2

    Returns list of tuples: [(sample_name, fragment_file), ...]
    Fragment file = prefix + '_atac_fragments.tsv.gz'
    """
    samples = []
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
            frag_file = f"{prefix}_atac_fragments.tsv.gz"
            samples.append((sample_name, frag_file))
    return samples

def main():
    parser = argparse.ArgumentParser(
        description="Consensus Peak Calling Script using MACS2 with samples_map"
    )
    parser.add_argument(
        "--samples_map", required=True,
        help="Sample map file. Format: '<prefix>: <sample_name>' per line"
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Output directory for BED and MACS2 peaks"
    )
    parser.add_argument(
        "-g", "--genome", default="mm",
        help="Genome size for MACS2 (default: mm)"
    )

    args = parser.parse_args()

    out_dir = args.output
    genome = args.genome

    print("=== STARTING PEAK CALLING SCRIPT ===", flush=True)
    print(f"Creating output directory: {out_dir}", flush=True)
    os.makedirs(out_dir, exist_ok=True)

    print(f"Reading samples map from: {args.samples_map}", flush=True)
    samples = read_samples_map(args.samples_map)
    print(f"Loaded {len(samples)} samples from map.", flush=True)

    for i, (sample, frag_path) in enumerate(samples):
        print(f">>> [{i+1}/{len(samples)}] Processing sample '{sample}' with fragments file '{frag_path}'", flush=True)

        if not os.path.exists(frag_path):
            print(f"  ERROR: Fragment file not found: {frag_path}", flush=True)
            continue

        bed_file = os.path.join(out_dir, f"{sample}.bed")

        # Convert fragments to BED WITHOUT SORTING (MACS2 can handle unsorted)
        print("  Step 1: Converting fragments to BED (no sort)...", flush=True)
        cmd_bed = f"zcat {frag_path} | awk 'BEGIN{{OFS=\"\\t\"}} {{print $1,$2,$3}}' > {bed_file}"
        print(f"  Running: {cmd_bed}", flush=True)

        try:
            subprocess.run(cmd_bed, shell=True, check=True, timeout=3600)  # 1 hour timeout
            file_size = os.path.getsize(bed_file) / (1024 * 1024)  # Size in MB
            print(f"  BED created: {file_size:.2f} MB", flush=True)
        except subprocess.TimeoutExpired:
            print("  TIMEOUT - skipping to next sample", flush=True)
            continue
        except Exception as e:
            print(f"  Conversion failed: {e}", flush=True)
            continue

        # Call MACS2
        print("  Step 2: Calling peaks with MACS2...", flush=True)
        cmd_macs = [
            "macs2", "callpeak",
            "-t", bed_file,
            "-f", "BED",
            "-n", sample,
            "--outdir", out_dir,
            "--nomodel",
            "--shift", "-100",
            "--extsize", "200",
            "-g", genome,
            "--keep-dup", "all"  # Important for ATAC-seq
        ]
        print(f"  Running: {' '.join(cmd_macs)}", flush=True)

        try:
            subprocess.run(cmd_macs, check=True, timeout=7200)  # 2 hour timeout
            print(f"  MACS2 completed for {sample}", flush=True)
        except subprocess.TimeoutExpired:
            print("  MACS2 TIMEOUT", flush=True)
        except Exception as e:
            print(f"  MACS2 failed: {e}", flush=True)

        print(f">>> Done {sample}", flush=True)
        print("-" * 50, flush=True)

    print("=== ALL SAMPLES PROCESSED ===", flush=True)


if __name__ == "__main__":
    main()

