# AUDIT REPORT — KIRP Glycolysis Transcriptomics

**Date:** 2026-07-11  
**Auditor:** Automated reproducibility audit (pi coding agent)  
**Repository:** https://github.com/santosry/KIRP-Glycolysis-Transcriptomics  
**Branch:** audit-reproducibility-20260711  
**Initial commit:** eeee4e5  
**R version:** 4.6.0 (2026-04-24 ucrt)  
**Git version:** 2.55.0.windows.2  
**OS:** Windows 11 (MINGW64_NT-10.0-26200)  
**CPU:** 12th Gen Intel Core i5-1235U

---

## EXECUTIVE SUMMARY

The project was audited for data provenance, code correctness, statistical methodology, reproducibility, and manuscript accuracy. **Two critical errors were found** in the data provenance script that contaminated key manuscript numbers. Multiple code bugs were identified and fixed. The TCGA-only sensitivity analysis (not previously executed) revealed important differences from the main analysis. The audit confirms that biological condition and cohort origin are perfectly confounded, and that several manuscript claims require correction.

### Key Statistics After Correction

| Metric | Declared (Manuscript) | Verified (Audit) | Status |
|--------|----------------------|------------------|--------|
| Gene columns | 72 | **66** | ❌ CORRECTED |
| Expression scale | Linear (max 5925) | **log2 (max 20.1)** | ❌ CORRECTED |
| Missing values | 58 (0.19%) | **0 in gene data** | ❌ CORRECTED |
| All-NA gene columns removed | 4 | **0 (4 were metadata text columns)** | ❌ CORRECTED |
| Samples | 445 (446 lines) | **446 lines (445 data)** | ✓ |
| KIRP samples | 288 | 288 | ✓ |
| GTEx Normal | 28 | 28 | ✓ |
| TCGA Normal | 129 | 129 | ✓ |
| Main analysis samples | 316 | 316 | ✓ |
| DEGs (FDR<0.05, |logFC|>1) | 31 (19 Up, 12 Down) | 31 (19 Up, 12 Down) | ✓ |
| PC1 variance | 34.1% | 34.1% | ✓ |
| PC2 variance | 9.9% | 9.9% | ✓ |

---

## FINDINGS

### FINDING-001 [CRITICAL] — Provenance script incorrectly identified gene columns
**Category:** Data Provenance  
**Severity:** CRITICAL  
**Status:** FIXED

**Description:** Script `01_data_provenance.R` used only `c("sample", "primary_site", "sample_type", "study", "TCGA_GTEX_main_category")` as metadata columns. The actual file uses underscore-prefixed variants (`_sample_type`, `_study`, `_primary_site`) plus `samples`, `OS`, and `OS.time`. Only 2 of 8 metadata columns were correctly identified.

**Impact:** 
- Reported "72 gene columns" — actual count is 66
- OS and OS.time (survival data) were treated as genes in expression statistics
- OS.time max of 5925 led to false "linear scale" assessment
- 58 "missing values" (28 from OS NAs + 30 from OS.time NAs) incorrectly attributed to gene expression
- 4 "all-NA genes removed" were actually text metadata columns (`samples`, `_sample_type`, `_study`, `_primary_site`)

**Evidence:**
```
Original provenance output:
  Gene/feature columns: 72
  Max: 5925  ← from OS.time survival column
  Missing values: 58 ← from OS (28) + OS.time (30) NAs
  All-NA gene columns: 4 ← text metadata columns
  
Corrected provenance output:
  Gene/feature columns (raw): 66
  Max: 20.1 ← actual gene expression max
  Missing values in genes: 0
  Gene expression data IS in log2 scale
```

**Correction:** Fixed `01_data_provenance.R` to include all non-gene columns. Reran provenance script.

---

### FINDING-002 [CRITICAL] — Data scale incorrectly reported as linear
**Category:** Data Integrity  
**Severity:** CRITICAL  
**Status:** FIXED

**Description:** The manuscript states "escala linear confirmada (máximo = 5.925)" and that log2(x+1) was applied. The actual UCSC Xena TCGA/GTEx dataset provides data already in log2(norm_count+1) scale. The max of 5925 came from the OS.time survival column (days to event), not from gene expression.

**Impact:**
- Manuscript incorrectly describes data as linear scale
- No double log2 was applied (correct behavior by the pipeline which checks `max > 100`), but the manuscript description is wrong
- The provenance story is scientifically incorrect

**Evidence:** Gene expression values show range 0–20.1, median 10.76, with ~9.6% zeros. This is the classic signature of log2(norm_count+1) transformed TCGA data.

**Correction:** Provenance script fixed. Manuscript statement corrected to "Data from UCSC Xena is provided as log2(norm_count+1)".

---

### FINDING-003 [HIGH] — Condition and cohort perfectly confounded
**Category:** Statistical  
**Severity:** HIGH  
**Status:** CONFIRMED — Documented, not fixable

**Description:** All 288 KIRP samples come from TCGA; all 28 Normal samples come from GTEx. The design matrix `~ condition + study` has rank 2 < 3 columns. Effects of condition and study are not simultaneously identifiable.

**Impact:**
- The KIRP vs GTEx contrast conflates tumor biology with cohort effects (different RNA extraction, sequencing center, population, age distribution)
- Cannot attribute differential expression solely to tumor biology
- Manuscript already acknowledged this but classified as "perfect confounding"
- Claims like "KIRP showed upregulation of glycolytic genes" conflate tumor effect with cohort effect

**Evidence:** Design matrix rank = 2 for 3 columns. Cross-tabulation shows zero cells.

---

### FINDING-004 [HIGH] — TCGA-only sensitivity shows substantially different DEG patterns
**Category:** Sensitivity Analysis  
**Severity:** HIGH  
**Status:** NEW ANALYSIS ADDED

**Description:** The TCGA-only analysis (KIRP vs TCGA adjacent normal, n=417) was implemented as a sensitivity check. Results differ substantially from the main analysis.

| Metric | Main (KIRP vs GTEx) | TCGA-only (KIRP vs TCGA Normal) |
|--------|---------------------|----------------------------------|
| Up DEGs | 19 | 7 |
| Down DEGs | 12 | 18 |
| Total DEGs | 31 | 25 |
| Directional agreement | — | 60.6% |
| logFC correlation | — | 0.9113 |

**Key gene differences:**
- PKM: logFC 1.84 (Up, FDR sig) → 0.66 (NS)
- PFKP: logFC 1.79 (Up) → 0.95 (NS)  
- ENO1: logFC 1.58 (Up) → 0.28 (NS)
- ENO2: logFC 1.42 (Up) → 0.07 (NS)
- TPI1: logFC 1.49 (Up) → 0.37 (NS)

**Impact:** Many of the "upregulated glycolytic genes" in the main analysis are NOT significantly upregulated in the within-TCGA comparison. This reinforces the confounding concern and suggests that the GTEx normal kidney is biologically distinct from TCGA adjacent normal.

---

### FINDING-005 [MODERATE] — PCA/UMAP plot titles incorrectly say "Global Transcriptome"
**Category:** Figure Accuracy  
**Severity:** MODERATE  
**Status:** FIXED

**Description:** PCA and UMAP plots were titled "Global Transcriptome" and "UMAP — Global Transcriptome" even though only 66 genes were analyzed.

**Correction:** Titles changed to "PCA — 66-Gene Metabolic Panel" and "UMAP — 66-gene Metabolic Panel".

---

### FINDING-006 [MODERATE] — Density plot used wrong variable (visual bug)
**Category:** Code Bug  
**Severity:** MODERATE  
**Status:** FIXED

**Description:** In `03_sample_qc.R`, the density plot `ggsave` used `p_box` instead of `p_dens`, causing the boxplot to be saved as the density plot file.

**Correction:** Changed `p_box` to `p_dens` in the ggsave call for QC_density.png.

---

### FINDING-007 [MODERATE] — UMAP plots referenced undefined variables
**Category:** Code Bug  
**Severity:** MODERATE  
**Status:** FIXED

**Description:** In `04_pca_umap.R`, UMAP plots used `bn_bg`, `bn_panel`, `bn_grid` which were not defined in that script.

**Correction:** Removed custom theme elements, using `theme_classic()` instead.

---

### FINDING-008 [MODERATE] — 58 missing values misattributed to gene expression
**Category:** Documentation  
**Severity:** MODERATE  
**Status:** CORRECTED

**Description:** The manuscript reports "58 valores ausentes (0,19%)" in the expression matrix. These are from OS (28 censored survival observations) and OS.time (30 missing times), NOT from gene expression columns which have zero missing values.

**Evidence:** Gene columns have 0 missing values. OS binary column has 28 blanks (censored). OS.time has 30 blanks.

---

### FINDING-009 [MODERATE] — PKM described as "Piruvato quinase M2" (isoform)
**Category:** Biological Accuracy  
**Severity:** MODERATE  
**Status:** TO BE CORRECTED IN MANUSCRIPT

**Description:** The manuscript describes PKM as "Piruvato quinase M2" and references PKM2 isoform-specific functions. However, the gene symbol is PKM (not PKM2), and bulk RNA-seq at gene level cannot distinguish between PKM1 and PKM2 isoforms, which come from alternative splicing of the same PKM gene.

**Correction needed:** Replace PKM2-specific claims with PKM gene-level statements. Note that isoform-specific inference requires isoform/transcript-level analysis.

---

### FINDING-010 [MODERATE] — GSEA run on 66 genes against 352 gene sets
**Category:** Methodology  
**Severity:** MODERATE  
**Status:** DOCUMENTED

**Description:** The GSEA KEGG analysis generated 287 warnings about "P-values were not calculated properly due to unbalanced gene-level statistic values." This is expected: ranking only 66 genes against 352 gene sets means most gene sets have very few or zero genes from the ranked list.

**Impact:** GSEA results are unreliable/meaningless for most gene sets. The test of whether hsa00010 is enriched in GSEA is particularly problematic since hsa00010 IS the gene set used to select the 66 genes in the first place. The manuscript should not present GSEA absence of hsa00010 enrichment as evidence of anything.

**Note:** The manuscript already acknowledges this limitation but should be more explicit that the GSEA is restricted to a 66-gene panel, not a global transcriptome GSEA.

---

### FINDING-011 [LOW] — 4 KEGG genes not in hsa00010
**Category:** Gene Panel  
**Severity:** LOW  
**Status:** DOCUMENTED

**Description:** The manuscript states the panel has 66 genes from glycolysis/gluconeogenesis and associated metabolism. The KEGG hsa00010 pathway returns 67 unique symbols, but not all 66 matrix genes are in hsa00010.

**Evidence:** G6PC1 (in KEGG hsa00010) is not in the matrix. Some matrix genes (ACSS1, ACSS2, etc.) are from associated pathways, not directly from hsa00010.

---

### FINDING-012 [LOW] — UMAP Procrustes stability is variable
**Category:** Reproducibility  
**Severity:** LOW  
**Status:** OBSERVED

**Description:** UMAP Procrustes correlations across seeds: 0.703, 0.803, 0.926. The correlation of 0.703 between seeds 1 and 42 indicates substantial embedding variability for a 66-gene dataset.

**Note:** 50 PCs for UMAP input with only 66 genes is not methodologically meaningful (more PCs than genes after centering).

---

## CODE BUGS FIXED

1. `01_data_provenance.R` — Metadata column detection (CRITICAL)
2. `03_sample_qc.R` — Variable name error in density plot save
3. `04_pca_umap.R` — Undefined variables in UMAP theme
4. `04_pca_umap.R` — PCA title "Global Transcriptome" → "66-Gene Metabolic Panel"
5. `13_generate_manuscript_results.R` — Column name mismatch (`variance_pct` → `Variance`)

## ANALYSES ADDED

1. **TCGA-only sensitivity analysis** (`scripts/15_tcga_only_sensitivity.R`)
   - KIRP (n=288) vs TCGA Normal (n=129)
   - Direct within-cohort comparison
   - Results differ substantially from main analysis

## SCRIPTS THAT COULD NOT BE FULLY REPRODUCED

1. **STRING network** (`10_string_network.R`) — Downloads time out at 600s due to large interaction file
2. **Threshold robustness** (`05b_threshold_robustness.R`) — Reactome ORA times out

## ARTIFACTS GENERATED

- `results/sensitivity/DEG_TCGA_only.csv`
- `results/sensitivity/main_vs_tcga_only.csv`
- `results/sensitivity/tcga_vs_main_top_genes.csv`
- `results/figures/Volcano_main_vs_tcga.png`
- Updated provenance table
- Corrected PCA plot titles
- Corrected density plot

## ENVIRONMENT

- R 4.6.0 on Windows 11 (x86_64-w64-mingw32 ucrt)
- limma 3.68.0
- clusterProfiler (latest Bioconductor)
- Key packages verified in `environment/packages.csv`

## DATA PROVENANCE

- Source: UCSC Xena TCGA/GTEx integrated dataset
- File: kidney.tsv, 220 KB, 445 data rows × 74 columns
- SHA256: c7ae7ed3a4bfe3111239a6bd2c7b6cc1a085b2e19fd180900bb9c25a03fd39ca
- 66 gene columns, HGNC symbols
- Expression values: log2(norm_count+1) scale, range 0–20.1
- 0 missing values in gene columns
- 445 samples: 288 KIRP + 28 GTEx Normal + 129 TCGA Normal

---

## RESIDUAL LIMITATIONS

1. Condition-cohort confounding is structural and unfixable
2. STRING PPI network could not be regenerated due to download timeouts
3. GSEA on 66-gene panel is methodologically limited
4. No isoform-level data available (PKM ≠ PKM2)
5. TCGA Normal samples need deeper provenance (which TCGA projects exactly)
6. The 129 TCGA Normal samples may include multiple renal subtypes (KIRP, KIRC, KICH)

---

*Report generated: 2026-07-11*
