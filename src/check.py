#!/usr/bin/env python3
"""
Check fragment distribution and compute meaningful stats
"""

import gzip
from collections import Counter
import numpy as np

files = [
    "13005-TH2_atac_fragments.tsv.gz",
    "13784-TH1_atac_fragments.tsv.gz",
    "13784-TH2_atac_fragments.tsv.gz"
]

print("📊 FRAGMENT DISTRIBUTION ANALYSIS")
print("=" * 70)

for f in files:
    print(f"\n📁 Sample: {f}")
    
    try:
        with gzip.open(f, 'rt') as gz:
            barcode_counts = Counter()
            fragments_read = 0
            
            for line in gz:
                if line.startswith('#'):
                    continue
                    
                cols = line.strip().split('\t')
                if len(cols) >= 4:
                    barcode = cols[3]
                    barcode_counts[barcode] += 1
                    fragments_read += 1
            
            if not barcode_counts:
                print("  ❗ No barcodes found!")
                continue
            
            # Convert to list of counts
            counts = list(barcode_counts.values())
            total_cells = len(counts)
            total_fragments = sum(counts)
            
            print(f"  Total fragments: {total_fragments:,}")
            print(f"  Total barcodes: {total_cells:,}")
            print(f"  Mean fragments per barcode: {np.mean(counts):.2f}")
            print(f"  Median fragments per barcode: {np.median(counts):.2f}")
            print(f"  Min fragments per barcode: {np.min(counts)}")
            print(f"  Max fragments per barcode: {np.max(counts)}")
            
            # Percentiles
            percentiles = [10, 25, 50, 75, 90, 95, 99]
            print(f"\n  Percentiles of fragments per barcode:")
            for p in percentiles:
                val = np.percentile(counts, p)
                print(f"    {p}th percentile: {val:.1f}")
            
            # Cells meeting common thresholds
            print(f"\n  Cells meeting thresholds:")
            thresholds = [10, 50, 100, 500, 1000]
            for t in thresholds:
                cells_above = sum(1 for c in counts if c >= t)
                percentage = (cells_above / total_cells) * 100
                print(f"    ≥{t} fragments: {cells_above:,} cells ({percentage:.1f}%)")
            
            # Distribution summary
            print(f"\n  Distribution summary:")
            print(f"    Top 10% of cells have ≥{np.percentile(counts, 90):.1f} fragments")
            print(f"    Bottom 10% of cells have ≤{np.percentile(counts, 10):.1f} fragments")
            
    except Exception as e:
        print(f"  ❌ ERROR: {e}")

print("\n" + "=" * 70)
print("RECOMMENDATION:")
print("Based on median fragments per cell, use --min_fragments_per_cb at median value")
print("or slightly above the 25th percentile.")
