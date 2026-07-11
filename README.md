# Central Carbon Metabolism Transcriptomics in KIRP

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R 4.6.0](https://img.shields.io/badge/R-4.6.0-276DC3)](https://www.r-project.org/)
[![v3.0.0](https://img.shields.io/badge/version-3.0.0-green)]()

**Transcriptomic profile of central carbon metabolism genes in papillary renal cell carcinoma (KIRP): paired tumor-adjacent analysis, sensitivity to reference tissue choice, and methodological implications.**

## Scientific Question

What is the transcriptomic profile of central carbon metabolism genes (glycolysis, pentose phosphate pathway, citrate cycle) in KIRP when analyzed with paired tumor-adjacent normal tissue as the primary comparator, and how sensitive are the results to the choice of reference tissue?

## Pathways Analyzed

| KEGG ID | Pathway | Genes in DB | Genes in Matrix |
|---------|---------|-------------|-----------------|
| hsa00010 | Glycolysis / Gluconeogenesis | 67 | 64 |
| hsa00030 | Pentose Phosphate Pathway | 31 | 30 |
| hsa00020 | Citrate Cycle (TCA) | 30 | 29 |
| **Union (unique)** | **Central Carbon Metabolism** | **110** | **106** |

## Key Findings (v3.0.0)

- **35/106 genes (33.0%)** differentially expressed in paired analysis (32 KIRP tumor-adjacent pairs, |log2FC|>1, FDR<0.05)
- **ALDOB** (log2FC = -8.66, CI: [-9.96, -7.36]) and **HK2** (+3.34, CI: [2.51, 4.17]) show largest magnitudes
- Camera gene set test: TCA genes coordinately suppressed (FDR = 0.0012)
- 7 genes shared across multiple pathways among the 35 unique DEGs
- **Paired vs TCGA-adjacent**: Lin's CCC = 0.990, MAE = 0.19, 100% directional agreement
- **Paired vs GTEx**: Lin's CCC = 0.795, MAE = 0.94, systematic bias = -0.79, only 60.2% directional agreement
- GTEx as reference tissue introduces systematic overestimation and directional discordance

## Comparator Design

| Comparator | Design | N | Role |
|------------|--------|---|------|
| Paired KIRP | 32 tumor-adjacent pairs, ~ patient + condition | 64 | Primary |
| KIRP vs Adjacent Normal | Unpaired, 288 tumors vs 32 KIRP-adjacent | 320 | Secondary (same cohort) |
| KIRP vs All TCGA Normal | Unpaired, 288 tumors vs 129 TCGA normals | 417 | Exploratory (mixed projects) |
| KIRP vs GTEx Normal | Cross-cohort, structurally confounded | 316 | Sensitivity only |

## Data

- **Source:** UCSC Xena — `TcgaTargetGtex_RSEM_Hugo_norm_count`
- **Samples:** 445 (288 KIRP, 28 GTEx Normal, 32 KIRP-adjacent normal, 97 other TCGA normal)
- **Genes:** 31,633 after filtering (58,581 raw, 26,948 low expression removed)
- **Scale:** log2(norm_count+1), confirmed by audit (max 20.1, median 11.15)
- **Access date:** 2026-07-11

## Reproduction

```r
# Clone
git clone https://github.com/santosry/KIRP-Glycolysis-Transcriptomics.git
cd KIRP-Glycolysis-Transcriptomics

# Download full kidney transcriptome (Python)
python3 scripts/download_full_transcriptome.py

# Run v3 pipeline
Rscript scripts/pipeline_v3.R
```

## Repository Structure

```
├── scripts/
│   ├── pipeline_v3.R                    # Primary analysis (paired, camera, concordance)
│   └── download_full_transcriptome.py   # Downloads & filters Xena data
├── results/v3/                         # v3.0.0 outputs
│   ├── tables/                          # DEGs, supplementary S1, discordant genes
│   ├── figures/                         # Volcano, Bland-Altman, PCA, paired plots
│   └── sessionInfo.txt                  # Complete session info
├── tests/                              # Automated tests (30/30 passing)
└── docs/                               # Documentation
```

## Limitations

- 32 pairs limits statistical power for small effects
- Bulk RNA-seq does not resolve cell types, isoforms, or tumor purity
- RNA abundance does not equal protein, enzyme activity, or metabolic flux
- GTEx comparison is structurally confounded (all tumors TCGA, all controls GTEx)
- KIRP is molecularly heterogeneous; pooled analysis does not capture subtypes
- No external validation cohort

## License

MIT License. See [LICENSE](LICENSE).

## Citation

See [CITATION.cff](CITATION.cff).
