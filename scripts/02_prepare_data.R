# 02_prepare_data.R
# PREPARO DA MATRIZ + FLUXO DE AMOSTRAS + QC INICIAL

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
raw_file <- file.path(repo_root, "data", "raw", "kidney.tsv")

if (!file.exists(raw_file)) {
  stop("Missing data/raw/kidney.tsv")
}

dir.create(file.path(repo_root, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "data", "metadata"), recursive = TRUE, showWarnings = FALSE)

# ── Import ──
kidney <- rio::import(raw_file)

# ── Sample flow tracking ──
sample_flow <- tibble(
  stage = character(),
  n_samples = integer(),
  description = character()
)

add_flow <- function(stage, n, desc) {
  sample_flow <<- rbind(sample_flow, tibble(stage = stage, n_samples = n, description = desc))
}

add_flow("01_imported", nrow(kidney), "Raw import from kidney.tsv")

# ── Remove duplicate samples ──
if ("sample" %in% colnames(kidney)) {
  dup_mask <- duplicated(kidney$sample)
  if (any(dup_mask)) {
    n_dup <- sum(dup_mask)
    message("Removing ", n_dup, " duplicate sample(s)")
    kidney <- kidney[!dup_mask, ]
  }
}
add_flow("02_dedup", nrow(kidney), "After duplicate removal")

# ── Required columns ──
required_cols <- c("sample", "primary_site", "sample_type", "study", "TCGA_GTEX_main_category")
missing_cols <- setdiff(required_cols, colnames(kidney))
if (length(missing_cols) > 0L) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# ── Separate metadata and expression ──
meta <- kidney[, required_cols, drop = FALSE]
gene_cols <- setdiff(colnames(kidney), required_cols)
expr <- as.matrix(kidney[, gene_cols, drop = FALSE])
storage.mode(expr) <- "numeric"

# ── Remove genes with zero variance ──
gene_vars <- apply(expr, 2, var, na.rm = TRUE)
zero_var_genes <- which(gene_vars == 0)
if (length(zero_var_genes) > 0) {
  message("Removing ", length(zero_var_genes), " zero-variance genes")
  expr <- expr[, -zero_var_genes, drop = FALSE]
}

# ── Remove all-NA genes ──
na_gene_cols <- which(apply(expr, 2, function(x) all(is.na(x))))
if (length(na_gene_cols) > 0) {
  message("Removing ", length(na_gene_cols), " all-NA genes")
  expr <- expr[, -na_gene_cols, drop = FALSE]
}

# ── Classify condition ──
meta$condition <- ifelse(
  meta$TCGA_GTEX_main_category == "GTEX Kidney", "Normal",
  ifelse(meta$TCGA_GTEX_main_category == "TCGA Kidney Papillary Cell Carcinoma", "KIRP", NA)
)
meta$condition <- factor(meta$condition, levels = c("Normal", "KIRP"))

# ── Condition × Study table (audit confounding) ──
confound_table <- as.data.frame(table(meta$condition, meta$study))
colnames(confound_table) <- c("condition", "study", "n")
rio::export(confound_table, file.path(repo_root, "data", "processed", "condition_by_study.csv"))
message("\nCondition × Study:")
print(confound_table)

# Check confounding
if (all(meta$study == "GTEX"[1] | is.na(meta$study))) {
  message("\n⚠ CONDIÇÃO BIOLÓGICA E ORIGEM DA COORTE ESTÃO PERFEITAMENTE CONFUNDIDAS")
  message("  KIRP → TCGA only, Normal → GTEx only")
  message("  This is a cross-cohort comparison. Interpret with caution.")
}

# ── Remove unclassified samples ──
unclassified <- sum(is.na(meta$condition))
if (unclassified > 0) {
  message("Removing ", unclassified, " unclassified samples")
  keep <- !is.na(meta$condition)
  meta <- meta[keep, , drop = FALSE]
  expr <- expr[keep, , drop = FALSE]
}
add_flow("03_classified", nrow(meta), "KIRP or Normal classified")

# ── Sample counts ──
n_kirp <- sum(meta$condition == "KIRP")
n_normal <- sum(meta$condition == "Normal")
message("\nFinal sample counts:")
message("  KIRP:  ", n_kirp)
message("  Normal: ", n_normal)
message("  Total:  ", nrow(meta))

add_flow("04_final", nrow(meta), paste0("Final: ", n_kirp, " KIRP + ", n_normal, " Normal"))

# ── Save ──
saveRDS(meta, file.path(repo_root, "data", "processed", "metadata.rds"))
saveRDS(expr, file.path(repo_root, "data", "processed", "expression_matrix.rds"))
rio::export(sample_flow, file.path(repo_root, "data", "processed", "sample_flow.csv"))

message("\n✓ Saved metadata.rds, expression_matrix.rds, sample_flow.csv")
