# Central Carbon Metabolism Transcriptomics in KIRP

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R 4.6.0](https://img.shields.io/badge/R-4.6.0-276DC3)](https://www.r-project.org/)
[![v2.0.0](https://img.shields.io/badge/version-2.0.0-green)]()

**Transcriptomic alterations in central carbon metabolism pathways in papillary renal cell carcinoma (KIRP): integrated analysis with multiple comparators and implications for precision health.**

---

## Scientific Question

What transcriptomic alterations are observed in interconnected central carbon metabolism pathways (glycolysis/gluconeogenesis, pentose phosphate pathway, citrate cycle) in KIRP, which remain robust across different comparator strategies, and what are the implications for interpreting omics evidence in precision health?

## Pathways Analyzed

| KEGG ID | Pathway | Genes in DB | Genes in Matrix |
|---------|---------|-------------|-----------------|
| hsa00010 | Glycolysis / Gluconeogenesis | 67 | 66 |
| hsa00030 | Pentose Phosphate Pathway | 31 | 11 |
| hsa00020 | Citrate Cycle (TCA) | 30 | 7 |
| **Union** | **Central Carbon Metabolism** | **110** | **66** |

## Key Findings (v2.0.0)

- **17/66 genes (25.8%)** are ROBUST — significant across all three comparators with consistent direction
- **ALDOB** (mean logFC = -8.11) and **HK2** (+3.30) are the most robustly altered
- ADH family genes (ADH1A, ADH1B, ADH1C, ADH4, ADH6) are consistently downregulated
- **PKM, PFKP, ENO1**: comparator-sensitive (significant only vs GTEx, not in paired/TCGA analyses)
- Paired analysis (32 pairs, gold standard) yields 7 Up / 19 Down, contrasting with GTEx: 19 Up / 12 Down
- Comparator choice substantially affects the gene list reported as differentially expressed

## Comparator Hierarchy

| Comparator | Design | N | Role |
|------------|--------|---|------|
| Paired KIRP | 32 tumor-normal pairs | 64 | Gold standard |
| KIRP vs TCGA Normal | Unpaired, same cohort | 417 | Primary unpaired |
| KIRP vs GTEx Normal | Cross-cohort, exploratory | 316 | Sensitivity |

## Data

- **Source:** UCSC Xena — TCGA/GTEx integrated dataset
- **Samples:** 445 (288 KIRP, 28 GTEx Normal, 32 KIRP-matched normal, 97 other TCGA normal)
- **Genes:** 66 (log2(norm_count+1) scale, max 20.1, 0 missing values)
- **SHA256:** `c7ae7ed3a4bfe3111239a6bd2c7b6cc1a085b2e19fd180900bb9c25a03fd39ca`

## Reproduction

```r
# Clone and run
git clone https://github.com/santosry/KIRP-Glycolysis-Transcriptomics.git
cd KIRP-Glycolysis-Transcriptomics

# Place kidney.tsv in data/raw/

# Run analysis
source("scripts/04_v2_main_analysis.R")
source("scripts/05_v2_visualizations.R")
```

## Repository Structure

```
├── data/raw/kidney.tsv          # Expression matrix
├── scripts/
│   ├── 01_column_reconciliation.R   # Resolves 72 vs 66 column discrepancy
│   ├── 02_tcga_normal_provenance.R  # Characterizes 129 TCGA Normal samples
│   ├── 03_kegg_pathway_genes.R      # KEGG pathway gene retrieval
│   ├── 04_v2_main_analysis.R        # Comparator hierarchy + robustness
│   └── 05_v2_visualizations.R       # Figures and tables
├── results/v2/                     # v2.0.0 outputs
├── docs/                           # Documentation
├── tests/                          # Automated tests
├── reports/                        # Audit reports
└── checksums/                      # SHA256SUMS
```

## Limitations

- PPP (11/31 genes) and TCA (7/30 genes) coverage is limited by the pre-filtered matrix
- RNA expression does not measure metabolic flux, enzyme activity, or protein levels
- Bulk RNA-seq does not resolve cell-type heterogeneity or isoforms (PKM1 vs PKM2)
- Cross-cohort (TCGA vs GTEx) comparison is structurally confounded
- No external validation cohort

## Precision Health Contribution

This study demonstrates that transcriptomic findings depend on comparator choice — a critical methodological consideration for interpreting omics evidence in precision health. Only 25.8% of genes were robust across comparators, underscoring the need for robustness assessment before biological or translational interpretation.

## License

MIT License. See [LICENSE](LICENSE).

## Citation

See [CITATION.cff](CITATION.cff) and [docs/MANUSCRIPT_VERSION_LOG.md](docs/MANUSCRIPT_VERSION_LOG.md).
