# Data Provenance

## Primary Data Source

**File:** `data/raw/kidney.tsv`  
**Name:** TCGA/GTEx integrated gene expression matrix  
**Source platform:** UCSC Xena (https://xenabrowser.net/)  
**Dataset:** GDC TCGA Kidney Papillary Cell Carcinoma (KIRP) combined with GTEx Kidney  
**Access date:** 2026-07-08  
**SHA256:** `c7ae7ed3a4bfe3111239a6bd2c7b6cc1a085b2e19fd180900bb9c25a03fd39ca`

## Data Characteristics

| Property | Value |
|----------|-------|
| File size | 220 KB |
| Rows (samples) | 445 |
| Total columns | 74 |
| Metadata columns | 8 (sample, samples, TCGA_GTEX_main_category, _sample_type, _study, _primary_site, OS, OS.time) |
| Gene columns | 66 |
| Gene identifiers | HGNC symbols |
| Expression scale | **log2(norm_count + 1)** (pre-transformed by UCSC Xena) |
| Expression range | 0 to 20.1 |
| Missing values (genes) | 0 |
| % zero expression | 9.6% |

## Sample Composition

| Group | n | Study | Sample Type |
|-------|---|-------|-------------|
| KIRP | 288 | TCGA | Primary Tumor |
| GTEx Normal | 28 | GTEx | Normal Tissue |
| TCGA Normal | 129 | TCGA | Solid Tissue Normal |

## Processing Steps

1. Downloaded from UCSC Xena as TSV
2. Metadata columns separated from gene expression columns
3. Underscore-prefixed columns renamed (_sample_type → sample_type, etc.)
4. Expression matrix extracted (66 gene columns, numeric)
5. **No log2 transformation applied** — data is already log2(norm_count+1)
6. Samples classified into KIRP, Normal_GTEx, and TCGA_Normal groups

## Technology

- RNA-seq: HTSeq counts
- Pipeline: UCSC Xena Toil RNA-seq recompute pipeline
- Data harmonization: TCGA and GTEx data were processed through a common pipeline
- Expression units: log2(norm_count + 1)

## Important Caveats

1. **Condition-cohort confounding:** All KIRP samples are from TCGA; all GTEx normals are from GTEx. The biological condition and cohort origin are perfectly confounded and not simultaneously identifiable.

2. **TCGA Normal samples:** 129 samples with blank TCGA_GTEX_main_category are classified as TCGA_Normal. These originate from TCGA projects but may include multiple renal carcinoma subtypes (KIRP, KIRC, KICH).

3. **This is a targeted panel, not a global transcriptome:** Only 66 genes from glycolysis/gluconeogenesis and associated metabolism pathways are included.

4. **Gene-level, not isoform-level:** The data measures gene-level expression. Isoform-specific inferences (e.g., PKM2 vs PKM1) are not supported.
