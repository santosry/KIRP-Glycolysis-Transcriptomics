# Reproducibility

## Environment

- **R version:** 4.6.0 (2026-04-24 ucrt)
- **Platform:** x86_64-w64-mingw32/x64
- **Package snapshot:** `environment/packages.csv`
- **Session info:** `environment/sessionInfo.txt`

## Reproduction Steps

### 1. Clone the repository
```bash
git clone https://github.com/santosry/KIRP-Glycolysis-Transcriptomics.git
cd KIRP-Glycolysis-Transcriptomics
```

### 2. Install dependencies
```r
source("scripts/00_environment.R")
```
This script checks for required packages and reports any missing ones.

### 3. Obtain the data
Download the TCGA/GTEx expression matrix from UCSC Xena:
1. Visit https://xenabrowser.net/datapages/
2. Select "GDC TCGA Kidney Papillary Cell Carcinoma (KIRP)" + "GTEx Kidney"
3. Download gene expression RNAseq (HTSeq - Counts)
4. Place as `data/raw/kidney.tsv`

Verify integrity:
```bash
sha256sum data/raw/kidney.tsv
# Expected: c7ae7ed3a4bfe3111239a6bd2c7b6cc1a085b2e19fd180900bb9c25a03fd39ca
```

### 4. Run the pipeline
```r
source("scripts/run_pipeline.R")
```

Or run individual steps:
```r
source("scripts/01_data_provenance.R")    # Audit and provenance
source("scripts/02_prepare_data.R")       # Classification, QC data prep
source("scripts/03_sample_qc.R")          # Sample-level QC
source("scripts/04_pca_umap.R")           # PCA and UMAP
source("scripts/05_differential_expression_global.R")  # limma DEG
source("scripts/06_hsa00010_targeted_analysis.R")      # hsa00010 analysis
source("scripts/07_ora_kegg.R")           # KEGG ORA
source("scripts/08_ora_reactome.R")       # Reactome ORA
source("scripts/09_rank_based_enrichment.R")  # GSEA
source("scripts/10_string_network.R")     # STRING PPI (requires internet)
source("scripts/15_tcga_only_sensitivity.R")  # TCGA-only sensitivity
```

### 5. Run tests
```r
testthat::test_dir("tests/testthat")
```

## Deterministic Components

- **limma eBayes:** Deterministic (given same inputs, produces identical results)
- **PCA:** Deterministic
- **ORA:** Deterministic
- **STRING PPI layout:** Seeded with `set.seed(1)` for reproducibility
- **UMAP:** Non-deterministic between seeds; Procrustes analysis provided for stability assessment

## Files Not Included in Repository

- `data/raw/kidney.tsv.gz` — Too large for GitHub (~70 MB compressed). Download from UCSC Xena.
- Manuscript files (`.Rmd`, `.pdf`, `.docx`) — NOT published per policy.

## Checksums

Output checksums are provided in `checksums/SHA256SUMS.txt` for verification of reproduced results.

## Known Limitations to Full Reproduction

1. **STRING network:** Requires internet connection and downloads ~500 MB of interaction data
2. **Network-dependent enrichment:** KEGG and STRING APIs may change over time
3. **UMAP:** Embedding varies between seeds; the pipeline records stability metrics

## CI Validation

GitHub Actions CI validates:
- Package availability
- All tests pass
- No manuscript files tracked
- No secrets exposed
