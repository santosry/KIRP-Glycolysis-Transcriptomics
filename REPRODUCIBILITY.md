# Reproducibility

## Environment

| Propriedade | Valor |
|-------------|-------|
| **R version** | 4.6.0 (2026-04-24 ucrt) |
| **Platform** | x86_64-w64-mingw32/x64 (Windows 11) |
| **Bioconductor** | 3.23 |
| **Package snapshot** | `environment/packages.csv` |
| **Session info** | `environment/sessionInfo.txt` |
| **renv lockfile** | `renv.lock` |

## Quick Start

### 1. Clone and setup

```bash
git clone https://github.com/santosry/KIRP-Glycolysis-Transcriptomics.git
cd KIRP-Glycolysis-Transcriptomics
```

### 2. Install R dependencies

```r
# Option A: Use renv (recommended)
install.packages("renv")
renv::restore()

# Option B: Manual installation
source("scripts/00_environment.R")
```

### 3. Download data

```bash
# Download full transcriptome (~178 MB, requires Python 3)
python3 scripts/download_full_transcriptome.py
```

### 4. Run pipeline

```r
# Primary analysis (paired, 32 pairs, camera, concordance)
source("scripts/pipeline_v3.R")

# Sensitivity addendum
source("scripts/pipeline_v3_addendum.R")

# 3D visualizations
source("scripts/16_3d_visualizations.R")       # Volcano 3D
source("scripts/16b_ppi_3d_correlation.R")      # PPI 3D
```

### 5. Run tests

```r
testthat::test_dir("tests/testthat")
# Expected: 30/30 passing
```

## Pipeline Steps (v3)

| Step | Script | Description |
|------|--------|-------------|
| Data download | `download_full_transcriptome.py` | Downloads 58,581 genes × 445 kidney samples from UCSC Xena |
| Full analysis | `pipeline_v3.R` | Loads, filters to 31,633 genes; paired DE (32 pairs); camera; concordance metrics; volcano plots |
| Sensitivity | `pipeline_v3_addendum.R` | Unpaired-64 model; treatment sensitivity; camera with inter.gene.cor=NA |
| 3D Volcano | `16_3d_visualizations.R` | Interactive 3D volcano plots per pathway (plotly) |
| 3D PPI | `16b_ppi_3d_correlation.R` | Gene co-expression network in 3D (Pearson |r| > 0.6, Fruchterman-Reingold) |

## Deterministic Components

| Component | Reproducibility | Notes |
|-----------|:---------------:|-------|
| limma eBayes | ✅ Deterministic | Given same inputs, produces identical results |
| PCA (prcomp) | ✅ Deterministic | |
| Pearson correlation | ✅ Deterministic | |
| Fisher's exact test (ORA) | ✅ Deterministic | |
| Camera gene set test | ✅ Deterministic | |
| Fruchterman-Reingold 3D | ✅ Deterministic | seed = 42 |
| KEGG annotations | ⚠️ Release-dependent | Frozen in `pathway_gene_membership.csv` |
| STRING PPI | ⚠️ Internet-dependent | API may change; alternative: correlation-based network |

## Files NOT Included in Repository

| File | Reason | How to obtain |
|------|--------|---------------|
| `data/raw/kidney_transcriptome.tsv` (~178 MB) | Too large for GitHub | Run `scripts/download_full_transcriptome.py` |
| `data/raw/*.gz` | Compressed duplicates | Run download scripts |
| `manuscrito_*` | Manuscript policy | Not published |
| `renv/library/` | Environment | Restore with `renv::restore()` |

## Checksums

Output checksums for all v3 results are in `results/v3/checksums_sha256.txt`. Verify with:

```bash
cd results/v3
sha256sum -c checksums_sha256.txt
```

## Known Limitations to Full Reproduction

1. **Transcriptome download:** `download_full_transcriptome.py` streams a ~2 GB gzip from UCSC Xena S3. Requires stable internet, 10-20 minutes.
2. **R package versions:** The `renv.lock` pins exact versions. If packages are removed from CRAN/Bioconductor, reproduction may fail.
3. **KEGG release:** Annotations were frozen from KEGG release 119.0 (July 2026). Future releases may alter pathway gene sets.
4. **Plotly HTML outputs:** The 3D HTML files require plotly.js 2.25.2 (bundled in `*_files/` directories). These must be kept alongside the HTML files.

## CI Validation

GitHub Actions CI (`.github/workflows/ci.yml`) validates:
- Package availability (all required packages installable)
- All 30 tests pass
- No manuscript files tracked in Git
- No secrets exposed

## Version History

| Version | Date | Key Changes |
|---------|------|-------------|
| v1.0.0 | 2026-07-08 | Initial: glycolysis only (hsa00010, 66 genes), KIRP vs GTEx |
| v2.0.0 | 2026-07-10 | Three pathways (106 genes), TCGA-only sensitivity, ORA+GSEA |
| v3.0.0 | 2026-07-11 | Full transcriptome (31,633 genes), paired primary analysis (32 pairs), camera, proper concordance (CCC, Bland-Altman, MAE) |
| v3.1.0 | 2026-07-12 | 3D interactive HTML visualizations, comprehensive README, AI agent instructions, ABNT references with DOI/access date |
