#!/usr/bin/env python3
"""
Download FULL kidney transcriptome from UCSC Xena (ALL genes, kidney samples only).
Streams TCGA/GTEx RSEM norm_count gzip, keeps all genes for kidney samples.
"""
import csv
import gzip
import os
import sys
import urllib.request
from io import TextIOWrapper

DATA_URL = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz"
KIDNEY_TSV = "data/raw/kidney.tsv"  # has sample IDs
OUTPUT_TSV = "data/raw/kidney_transcriptome.tsv"

print("=== Step 1: Reading kidney sample IDs ===")
kidney_samples = set()
with open(KIDNEY_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        kidney_samples.add(row["sample"])
print(f"Kidney samples: {len(kidney_samples)}")

print("=== Step 2: Downloading & filtering full matrix (~2 GB gzip) ===")
req = urllib.request.urlopen(DATA_URL)
decompressor = gzip.GzipFile(fileobj=req)
reader = TextIOWrapper(decompressor, encoding="utf-8", errors="replace")

header_line = reader.readline().strip()
all_columns = header_line.split("\t")
print(f"Total columns: {len(all_columns)}")

# Find kidney sample column indices
sample_indices = []
sample_names = []
for i, col_name in enumerate(all_columns):
    if col_name in kidney_samples:
        sample_indices.append(i)
        sample_names.append(col_name)

print(f"Kidney columns found: {len(sample_indices)}")

if len(sample_indices) == 0:
    print("ERROR: No kidney samples found!")
    sys.exit(1)

print("=== Step 3: Streaming all genes (this will take 10-20 min) ===")
# We'll write in batches to manage memory
os.makedirs(os.path.dirname(OUTPUT_TSV), exist_ok=True)
fout = open(OUTPUT_TSV, "w", buffering=1024*1024)
writer = csv.writer(fout, delimiter="\t", lineterminator="\n")

# Write header: sample + gene symbols
writer.writerow(["gene"] + sample_names)

total_lines = 0
kept_lines = 0
for line in reader:
    total_lines += 1
    if total_lines % 5000 == 0:
        print(f"  Processed {total_lines} genes, kept {kept_lines}...")
    
    parts = line.strip().split("\t")
    if len(parts) < len(all_columns):
        continue
    
    gene_symbol = parts[0]
    # Skip genes with empty symbols or problematic names
    if not gene_symbol or gene_symbol == "?" or "?" in gene_symbol:
        continue
    
    # Extract expression values for kidney samples
    values = [gene_symbol] + [parts[i] if i < len(parts) else "" for i in sample_indices]
    writer.writerow(values)
    kept_lines += 1

fout.close()
reader.close()
req.close()

print(f"\n=== Complete ===")
print(f"Genes processed: {total_lines}")
print(f"Genes kept: {kept_lines}")

# Verify
with open(OUTPUT_TSV) as f:
    header = f.readline().strip().split("\t")
    nsamples = len(header) - 1
print(f"Output: {nsamples} samples x {kept_lines} genes")

import os
size_mb = os.path.getsize(OUTPUT_TSV) / 1e6
print(f"File size: {size_mb:.0f} MB")
