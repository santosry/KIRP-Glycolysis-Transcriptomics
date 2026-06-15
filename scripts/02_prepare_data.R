suppressPackageStartupMessages({
  library(dplyr)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
raw_file <- file.path(repo_root, "data", "raw", "kidney.tsv")

if (!file.exists(raw_file)) {
  stop("Missing data/raw/kidney.tsv. Place the UCSC Xena TCGA/GTEx expression matrix there before running the pipeline.")
}

dir.create(file.path(repo_root, "data", "processed"), recursive = TRUE, showWarnings = FALSE)

kidney <- rio::import(raw_file)
stopifnot(is.data.frame(kidney))

if ("samples" %in% colnames(kidney)) kidney$samples <- NULL
names(kidney) <- sub("^_", "", names(kidney))

required_cols <- c("sample", "primary_site", "sample_type", "study", "TCGA_GTEX_main_category")
missing_cols <- setdiff(required_cols, colnames(kidney))
if (length(missing_cols) > 0L) {
  stop("Missing required columns in kidney.tsv: ", paste(missing_cols, collapse = ", "))
}

required_groups <- c("GTEX Kidney", "TCGA Kidney Papillary Cell Carcinoma")
missing_groups <- setdiff(required_groups, unique(kidney$TCGA_GTEX_main_category))
if (length(missing_groups) > 0L) {
  stop("Missing required biological groups: ", paste(missing_groups, collapse = ", "))
}

meta <- kidney[, required_cols, drop = FALSE]
expr <- kidney[, setdiff(colnames(kidney), required_cols), drop = FALSE]
expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

meta$condition <- ifelse(
  meta$TCGA_GTEX_main_category == "GTEX Kidney", "Normal",
  ifelse(meta$TCGA_GTEX_main_category == "TCGA Kidney Papillary Cell Carcinoma", "KIRP", NA)
)
meta$condition <- factor(meta$condition, levels = c("Normal", "KIRP"))

keep_samples <- !is.na(meta$condition)
meta <- meta[keep_samples, , drop = FALSE]
expr <- expr[keep_samples, , drop = FALSE]

if (all(expr >= 0, na.rm = TRUE) && all(abs(expr - round(expr)) < .Machine$double.eps^0.5, na.rm = TRUE)) {
  warning("Expression matrix looks integer-like. Confirm it is log2(norm_count + 1).")
}

saveRDS(meta, file.path(repo_root, "data", "processed", "metadata.rds"))
saveRDS(expr, file.path(repo_root, "data", "processed", "expression_matrix.rds"))

qc_tables <- list(
  condition_by_sample_type = as.data.frame(table(meta$condition, meta$sample_type)),
  condition_by_study = as.data.frame(table(meta$condition, meta$study))
)
saveRDS(qc_tables, file.path(repo_root, "data", "processed", "qc_tables.rds"))

message("Saved processed metadata and expression matrix in data/processed/")
