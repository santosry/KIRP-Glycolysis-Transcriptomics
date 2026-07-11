#!/usr/bin/env python3
"""
Download full expression matrix for all 110 central carbon metabolism genes
from UCSC Xena TCGA/GTEx dataset. Streams the gzip file and filters on the fly.
"""
import csv
import gzip
import os
import sys
import urllib.request
from io import TextIOWrapper

# ── Configuration ──
DATA_URL = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz"
GENE_MEMBERSHIP = "results/tables/pathway_gene_membership.csv"
KIDNEY_TSV = "data/raw/kidney.tsv"
OUTPUT_TSV = "data/raw/kidney_full.tsv"
OUTPUT_GZ = "data/raw/kidney_full.tsv.gz"

# ── Step 1: Get all pathway genes ──
print("=== Step 1: Reading pathway genes ===")
pathway_genes = set()
with open(GENE_MEMBERSHIP) as f:
    for row in csv.DictReader(f):
        pathway_genes.add(row["gene"])
print(f"Pathway genes: {len(pathway_genes)}")

# ── Step 2: Get sample IDs ──
print("=== Step 2: Reading sample IDs ===")
kidney_samples = set()
with open(KIDNEY_TSV) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        kidney_samples.add(row["sample"])
print(f"Kidney samples: {len(kidney_samples)}")

# ── Step 3: Download and filter ──
print("=== Step 3: Downloading full matrix (this will take a while) ===")
print(f"URL: {DATA_URL}")

# The matrix is genes×samples (first col = gene symbol, rest = sample IDs)
req = urllib.request.urlopen(DATA_URL)
decompressor = gzip.GzipFile(fileobj=req)
reader = TextIOWrapper(decompressor, encoding="utf-8")

# Read header to get sample column indices
header_line = reader.readline().strip()
all_columns = header_line.split("\t")
print(f"Total columns in full matrix: {len(all_columns)}")

# Find which columns correspond to our kidney samples
sample_col_indices = []
sample_col_names = []
for i, col_name in enumerate(all_columns):
    if col_name in kidney_samples:
        sample_col_indices.append(i)
        sample_col_names.append(col_name)
    elif col_name == "sample" or col_name == "gene":
        pass  # skip metadata col

print(f"Kidney sample columns found: {len(sample_col_indices)}")

if len(sample_col_indices) == 0:
    print("ERROR: No kidney samples found in full matrix!")
    print("First 10 columns:", all_columns[:10])
    print("First 5 kidney samples:", list(kidney_samples)[:5])
    sys.exit(1)

# Find gene column (first column)
gene_col_idx = 0  # "sample" column in Xena format = gene symbols

# ── Step 4: Stream and filter ──
print("=== Step 4: Filtering genes (streaming) ===")
output_rows = []
genes_found = 0
total_lines = 0

for line in reader:
    total_lines += 1
    if total_lines % 5000 == 0:
        print(f"  Processed {total_lines} genes, found {genes_found}...")
    
    parts = line.strip().split("\t")
    if len(parts) < len(all_columns):
        continue
    
    gene = parts[gene_col_idx]
    if gene in pathway_genes:
        # Extract expression values for kidney samples
        values = [parts[i] for i in sample_col_indices]
        output_rows.append([gene] + values)
        genes_found += 1
        
        if genes_found == len(pathway_genes):
            print(f"  All {len(pathway_genes)} genes found! Stopping early.")
            break

reader.close()
req.close()

print(f"Genes found: {genes_found} / {len(pathway_genes)}")
print(f"Total lines processed: {total_lines}")

# ── Step 5: Write output ──
print("=== Step 5: Writing output ===")
# Transpose to samples×genes format (matching original kidney.tsv format)
# Input: genes × samples → Output: samples × genes
os.makedirs(os.path.dirname(OUTPUT_TSV), exist_ok=True)

with open(OUTPUT_TSV, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    # Header: sample, then gene names
    writer.writerow(["sample"] + [row[0] for row in output_rows])
    # Rows: one per sample
    for si, sample_name in enumerate(sample_col_names):
        row_data = [sample_name] + [output_rows[gi][si + 1] for gi in range(len(output_rows))]
        writer.writerow(row_data)

print(f"Written: {OUTPUT_TSV}")
print(f"Samples: {len(sample_col_names)}, Genes: {genes_found}")

# Also gzip it
import gzip as gz
with open(OUTPUT_TSV, "rb") as fin:
    with gz.open(OUTPUT_GZ, "wb") as fout:
        fout.write(fin.read())
print(f"Gzipped: {OUTPUT_GZ}")

# ── Verify ──
print("=== Verification ===")
with open(OUTPUT_TSV) as f:
    header = f.readline().strip().split("\t")
    n_genes = len(header) - 1
    n_samples = sum(1 for _ in f)
print(f"Output dimensions: {n_samples} samples × {n_genes} genes")

# Missing genes
missing = pathway_genes - set(header[1:])
if missing:
    print(f"Missing genes ({len(missing)}): {sorted(missing)}")
else:
    print("All genes found!")
