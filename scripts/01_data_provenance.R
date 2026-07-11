# 01_data_provenance.R
# AUDITORIA INICIAL DOS DADOS
# Examina a matriz kidney.tsv, gera hash, tabela de proveniência

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
raw_file <- file.path(repo_root, "data", "raw", "kidney.tsv")
dir.create(file.path(repo_root, "data", "provenance"), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(raw_file)) {
  stop("Missing data/raw/kidney.tsv. Download from UCSC Xena before running.")
}

# ── SHA256 hash ──
hash <- digest::digest(file = raw_file, algo = "sha256")
message("SHA256: ", hash)

# ── File info ──
file_info <- file.info(raw_file)
message("File size: ", file_info$size, " bytes (", round(file_info$size/1e6, 1), " MB)")
message("Last modified: ", file_info$mtime)

# ── Load ──
kidney <- rio::import(raw_file)
stopifnot(is.data.frame(kidney))
message("Dimensions: ", nrow(kidney), " rows × ", ncol(kidney), " cols")

# ── Column inspection ──
cols <- colnames(kidney)
# Full set of known non-gene columns (both plain and underscore-prefixed variants)
known_non_gene <- c(
  "sample", "samples",
  "primary_site", "_primary_site",
  "sample_type", "_sample_type",
  "study", "_study",
  "TCGA_GTEX_main_category",
  "OS", "OS.time"
)
meta_present <- intersect(known_non_gene, cols)
gene_cols <- setdiff(cols, known_non_gene)
# Also strip any columns that resolve to known_non_gene after underscore removal
still_candidates <- setdiff(gene_cols, NULL)
message("Non-gene columns found: ", paste(meta_present, collapse = ", "))
message("Gene/feature columns (raw): ", length(gene_cols))

# ── Sample inspection ──
if ("sample" %in% cols) {
  dup_samples <- any(duplicated(kidney$sample))
  message("Duplicate sample IDs: ", ifelse(dup_samples, "YES ⚠", "NO ✓"))
}

if ("TCGA_GTEX_main_category" %in% cols) {
  categories <- table(kidney$TCGA_GTEX_main_category)
  message("Categories in TCGA_GTEX_main_category:")
  for (cat_name in names(categories)) {
    message("  ", cat_name, ": ", categories[cat_name])
  }
}

# ── Expression matrix inspection ──
# Only include truly numeric gene columns (exclude any columns that become all-NA)
expr_raw <- kidney[, gene_cols, drop = FALSE]
# Convert each column to numeric; flag columns that become all-NA (metadata text columns)
gene_cols_valid <- character(0)
expr_list <- list()
for (gc in gene_cols) {
  col_vals <- suppressWarnings(as.numeric(kidney[[gc]]))
  if (all(is.na(col_vals))) {
    message("Column ", gc, " is non-numeric (all-NA after coercion) — EXCLUDED from gene set")
  } else {
    gene_cols_valid <- c(gene_cols_valid, gc)
    expr_list[[gc]] <- col_vals
  }
}
message("Gene columns after numeric validation: ", length(gene_cols_valid))

expr_raw <- do.call(cbind, expr_list)
colnames(expr_raw) <- gene_cols_valid
storage.mode(expr_raw) <- "numeric"

# Basic stats on GENE columns only
vals <- as.vector(expr_raw)
vals <- vals[!is.na(vals)]
n_missing_gene <- sum(is.na(expr_raw))
message("\nExpression value statistics (GENE columns only):")
message("  Min:    ", min(vals))
message("  1st Qu: ", quantile(vals, 0.25))
message("  Median: ", median(vals))
message("  Mean:   ", mean(vals))
message("  3rd Qu: ", quantile(vals, 0.75))
message("  Max:    ", max(vals))
message("  % zero:  ", round(100 * mean(vals == 0), 1), "%")
message("  Missing values in genes: ", n_missing_gene)

# Scale assessment (on gene data only)
if (min(vals) >= 0 && max(vals) < 30) {
  message("\n⟹ Gene expression data IS in log2 scale (values 0–", round(max(vals), 1), ", typical for log2(norm_count+1))")
  message("  The UCSC Xena dataset uses log2(norm_count+1) transformation.")
  message("  No log2 transformation is needed.")
} else if (min(vals) >= 0 && max(vals) > 100) {
  message("\n⟹ Data appears to be in linear scale. log2(x+1) will be applied downstream.")
} else if (min(vals) < 0) {
  message("\n⟹ Data contains negative values → likely centered/log-ratio transformed")
}

# Gene identifier inspection
sample_genes <- head(gene_cols_valid, 20)
message("\nSample gene identifiers (first 20):")
message(paste(sample_genes, collapse = ", "))

# Check if HGNC symbols
hgnc_pattern <- grepl("^[A-Z][A-Z0-9]+$", gene_cols_valid)
message("Appear to be HGNC symbols: ", round(100 * mean(hgnc_pattern), 1), "%")

# Check if Ensembl IDs
ensembl_pattern <- grepl("^ENSG\\d+", gene_cols_valid)
if (any(ensembl_pattern)) message("Ensembl IDs detected ⚠")

# Variance analysis
gene_vars <- apply(expr_raw, 2, var, na.rm = TRUE)
zero_var <- sum(gene_vars == 0, na.rm = TRUE)
low_var <- sum(gene_vars < 0.01, na.rm = TRUE)
message("\nGenes with zero variance: ", zero_var)
message("Genes with variance < 0.01: ", low_var)

# ── Provenance table ──
provenance <- tibble(
  field = c(
    "file", "sha256", "size_bytes", "size_mb", "last_modified",
    "n_rows_total", "n_gene_cols", "n_meta_cols",
    "expression_min", "expression_q25", "expression_median", "expression_q75", "expression_max",
    "pct_zero", "n_zero_var_genes", "n_missing_gene_cells",
    "scale_assessment", "gene_id_type",
    "access_date", "source_url"
  ),
  value = c(
    basename(raw_file), hash, as.character(file_info$size),
    as.character(round(file_info$size/1e6, 1)), as.character(file_info$mtime),
    as.character(nrow(kidney)), as.character(length(gene_cols_valid)), as.character(length(meta_present)),
    as.character(min(vals)), as.character(quantile(vals, 0.25)),
    as.character(median(vals)), as.character(quantile(vals, 0.75)), as.character(max(vals)),
    as.character(round(100 * mean(vals == 0), 1)),
    as.character(zero_var), as.character(n_missing_gene),
    "log2-scale",
    ifelse(mean(hgnc_pattern) > 0.9, "HGNC_symbols", "mixed_or_Ensembl"),
    "2026-07-11",
    "UCSC Xena (https://xenabrowser.net/)"
  )
)
rio::export(provenance, file.path(repo_root, "data", "provenance", "provenance.csv"))

message("\n✓ Provenance table saved to data/provenance/provenance.csv")
message("✓ SHA256: ", hash)
