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
meta_cols <- c("sample", "primary_site", "sample_type", "study", "TCGA_GTEX_main_category")
meta_present <- intersect(meta_cols, cols)
gene_cols <- setdiff(cols, meta_cols)
message("Metadata columns found: ", paste(meta_present, collapse = ", "))
message("Gene/feature columns: ", length(gene_cols))

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
expr_raw <- kidney[, gene_cols, drop = FALSE]
expr_raw <- as.matrix(expr_raw)
storage.mode(expr_raw) <- "numeric"

# Remove columns that are all NA
na_cols <- which(apply(expr_raw, 2, function(x) all(is.na(x))))
if (length(na_cols) > 0) {
  message("All-NA gene columns: ", length(na_cols), " (removing)")
  expr_raw <- expr_raw[, -na_cols, drop = FALSE]
  gene_cols <- gene_cols[-na_cols]
}

# Remove rows that are all NA
na_rows <- which(apply(expr_raw, 1, function(x) all(is.na(x))))
if (length(na_rows) > 0) {
  message("All-NA sample rows: ", length(na_rows))
}

# Basic stats
vals <- as.vector(expr_raw)
vals <- vals[!is.na(vals)]
message("\nExpression value statistics:")
message("  Min:    ", min(vals))
message("  1st Qu: ", quantile(vals, 0.25))
message("  Median: ", median(vals))
message("  Mean:   ", mean(vals))
message("  3rd Qu: ", quantile(vals, 0.75))
message("  Max:    ", max(vals))
message("  % integer-like: ", round(100 * mean(abs(vals - round(vals)) < 1e-10), 1), "%")

# Scale assessment
if (min(vals) >= 0 && max(vals) < 30) {
  message("\n⟹ Data appears to be in log2 scale (values 0–30, typical for log2(norm_count+1))")
} else if (min(vals) >= 0 && max(vals) > 100) {
  message("\n⟹ Data may be in linear scale (wide range, min ≥ 0)")
} else if (min(vals) < 0) {
  message("\n⟹ Data contains negative values → likely centered/log-ratio transformed")
}

# Gene identifier inspection
sample_genes <- head(gene_cols, 20)
message("\nSample gene identifiers (first 20):")
message(paste(sample_genes, collapse = ", "))

# Check if HGNC symbols
hgnc_pattern <- grepl("^[A-Z][A-Z0-9]+$", sample_genes)
message("Appear to be HGNC symbols: ", round(100 * mean(hgnc_pattern), 1), "%")

# Check if Ensembl IDs
ensembl_pattern <- grepl("^ENSG\\d+", sample_genes)
if (any(ensembl_pattern)) message("Ensembl IDs detected ⚠")

# Check for version suffixes (.digit)
has_version <- grepl("\\.\\d+$", sample_genes)
if (any(has_version)) message("Gene IDs with version suffix detected ⚠")

# Variance analysis
gene_vars <- apply(expr_raw, 2, var, na.rm = TRUE)
zero_var <- sum(gene_vars == 0, na.rm = TRUE)
low_var <- sum(gene_vars < 0.01, na.rm = TRUE)
message("\nGenes with zero variance: ", zero_var)
message("Genes with variance < 0.01: ", low_var)

# Missing values
missing_cells <- sum(is.na(expr_raw))
message("Missing values: ", missing_cells, " (", round(100*missing_cells/length(expr_raw), 3), "%)")

# Inf values
inf_cells <- sum(is.infinite(expr_raw))
message("Infinite values: ", inf_cells)

# ── Provenance table ──
provenance <- tibble(
  field = c(
    "file", "sha256", "size_bytes", "size_mb", "last_modified",
    "n_rows_total", "n_gene_cols", "n_meta_cols",
    "expression_min", "expression_q25", "expression_median", "expression_q75", "expression_max",
    "pct_integer_like", "n_zero_var_genes", "n_missing_cells", "n_inf_cells",
    "scale_assessment", "gene_id_type",
    "access_date", "source_url"
  ),
  value = c(
    basename(raw_file), hash, as.character(file_info$size),
    as.character(round(file_info$size/1e6, 1)), as.character(file_info$mtime),
    as.character(nrow(kidney)), as.character(length(gene_cols)), as.character(length(meta_present)),
    as.character(min(vals)), as.character(quantile(vals, 0.25)),
    as.character(median(vals)), as.character(quantile(vals, 0.75)), as.character(max(vals)),
    as.character(round(100 * mean(abs(vals - round(vals)) < 1e-10), 1)),
    as.character(zero_var), as.character(missing_cells), as.character(inf_cells),
    ifelse(min(vals) >= 0 && max(vals) < 30, "log2-scale", "unknown"),
    ifelse(mean(hgnc_pattern) > 0.9, "HGNC_symbols", "mixed_or_Ensembl"),
    "2026-07-08",
    "UCSC Xena (https://xenabrowser.net/)"
  )
)
rio::export(provenance, file.path(repo_root, "data", "provenance", "provenance.csv"))

message("\n✓ Provenance table saved to data/provenance/provenance.csv")
message("✓ SHA256: ", hash)
