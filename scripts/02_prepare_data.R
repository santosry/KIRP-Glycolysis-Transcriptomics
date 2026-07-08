# 02_prepare_data.R
# CLASSIFICAÇÃO DE AMOSTRAS + AUDITORIA DE CONFUNDIMENTO

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
raw_file <- file.path(repo_root, "data", "raw", "kidney.tsv")
if (!file.exists(raw_file)) stop("Missing data/raw/kidney.tsv. Run 00_download_data.R first or download manually.")

dir.create(file.path(repo_root, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "data", "metadata"), recursive = TRUE, showWarnings = FALSE)

kidney <- rio::import(raw_file)
message("Imported: ", nrow(kidney), " rows × ", ncol(kidney), " cols")

# ── Required columns ──
required_cols <- c("sample", "primary_site", "sample_type", "study", "TCGA_GTEX_main_category")
missing_cols <- setdiff(required_cols, colnames(kidney))
if (length(missing_cols) > 0L) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

# ── Deduplicate ──
if (any(duplicated(kidney$sample))) {
  n_dup <- sum(duplicated(kidney$sample))
  message("Removing ", n_dup, " duplicate samples")
  kidney <- kidney[!duplicated(kidney$sample), ]
}

# ── CLASSIFICATION AUDIT ──
# Export all cross-tabulations BEFORE classifying
audit_dir <- file.path(repo_root, "data", "processed")
meta_cols <- required_cols

# Study distribution
t_study <- as.data.frame(table(kidney$study))
rio::export(t_study, file.path(audit_dir, "audit_table_study.csv"))

# sample_type distribution
t_stype <- as.data.frame(table(kidney$sample_type))
rio::export(t_stype, file.path(audit_dir, "audit_table_sample_type.csv"))

# TCGA_GTEX_main_category
t_cat <- as.data.frame(table(kidney$TCGA_GTEX_main_category))
rio::export(t_cat, file.path(audit_dir, "audit_table_category.csv"))

# study × sample_type
t_ss <- as.data.frame(table(kidney$study, kidney$sample_type))
rio::export(t_ss, file.path(audit_dir, "audit_table_study_x_sample_type.csv"))

# category × sample_type
t_cs <- as.data.frame(table(kidney$TCGA_GTEX_main_category, kidney$sample_type))
rio::export(t_cs, file.path(audit_dir, "audit_table_category_x_sample_type.csv"))

# study × category × sample_type
t_scs <- as.data.frame(table(kidney$study, kidney$TCGA_GTEX_main_category, kidney$sample_type))
rio::export(t_scs, file.path(audit_dir, "audit_table_study_x_category_x_sample_type.csv"))

# ── Build classification audit table ──
class_audit <- kidney |>
  group_by(study, TCGA_GTEX_main_category, sample_type) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(
    proposed_condition = case_when(
      TCGA_GTEX_main_category == "GTEX Kidney" ~ "Normal_GTEx",
      TCGA_GTEX_main_category == "TCGA Kidney Papillary Cell Carcinoma" &
        grepl("Primary|Tumor", sample_type, ignore.case = TRUE) ~ "KIRP",
      TCGA_GTEX_main_category == "TCGA Kidney Papillary Cell Carcinoma" &
        grepl("Solid Tissue Normal|Normal", sample_type, ignore.case = TRUE) ~ "TCGA_Normal",
      TRUE ~ "UNCLASSIFIED"
    ),
    reason = case_when(
      TCGA_GTEX_main_category == "GTEX Kidney" ~ "GTEx normal kidney",
      TCGA_GTEX_main_category == "TCGA Kidney Papillary Cell Carcinoma" &
        grepl("Primary|Tumor", sample_type, ignore.case = TRUE) ~ "TCGA-KIRP primary tumor",
      TCGA_GTEX_main_category == "TCGA Kidney Papillary Cell Carcinoma" &
        grepl("Solid Tissue Normal|Normal", sample_type, ignore.case = TRUE) ~ "TCGA-KIRP adjacent normal",
      TRUE ~ "Unknown classification"
    )
  )

rio::export(class_audit, file.path(audit_dir, "sample_classification_audit.csv"))
message("\nClassification audit:")
print(as.data.frame(class_audit))

# ── Separate metadata and expression ──
meta <- kidney[, required_cols, drop = FALSE]
gene_cols <- setdiff(colnames(kidney), required_cols)
expr <- as.matrix(kidney[, gene_cols, drop = FALSE])
storage.mode(expr) <- "numeric"

# ── Build condition ──
meta <- meta |>
  left_join(class_audit |> select(study, TCGA_GTEX_main_category, sample_type, proposed_condition),
            by = c("study", "TCGA_GTEX_main_category", "sample_type"))

# Filter to classified samples
meta <- meta |> filter(proposed_condition != "UNCLASSIFIED")
expr <- expr[match(meta$sample, kidney$sample), , drop = FALSE]

n_kirp     <- sum(meta$proposed_condition == "KIRP")
n_gtx_norm <- sum(meta$proposed_condition == "Normal_GTEx")
n_tcga_norm <- sum(meta$proposed_condition == "TCGA_Normal")

message("\nSample counts after classification:")
message("  KIRP (TCGA primary tumor): ", n_kirp)
message("  Normal GTEx: ", n_gtx_norm)
message("  TCGA Normal (adjacent): ", n_tcga_norm)

# ── CONFOUNDING AUDIT ──
# For main analysis: KIRP vs Normal_GTEx
main_meta <- meta |> filter(proposed_condition %in% c("KIRP", "Normal_GTEx"))
main_meta$condition <- factor(main_meta$proposed_condition, levels = c("Normal_GTEx", "KIRP"))

# condition × study
confound_table <- as.data.frame(table(main_meta$condition, main_meta$study))
colnames(confound_table) <- c("condition", "study", "n")
rio::export(confound_table, file.path(audit_dir, "condition_by_study.csv"))

message("\nCondition × Study:")
print(confound_table)

# Check design matrix rank
design_check <- model.matrix(~ condition + study, data = main_meta)
design_rank <- qr(design_check)$rank
design_ncol <- ncol(design_check)

message("Design matrix rank (condition + study): ", design_rank, " / ", design_ncol, " cols")
if (design_rank < design_ncol) {
  message("⚠ Effects of condition and study are NOT simultaneously identifiable.")
}

# Perfect confounding check
studies_kirp <- unique(main_meta$study[main_meta$condition == "KIRP"])
studies_normal <- unique(main_meta$study[main_meta$condition == "Normal_GTEx"])
if (length(intersect(studies_kirp, studies_normal)) == 0) {
  message("\n⚠ CONDIÇÃO BIOLÓGICA E ORIGEM DA COORTE ESTÃO PERFEITAMENTE CONFUNDIDAS")
  message("  KIRP samples come from: ", paste(studies_kirp, collapse = ", "))
  message("  Normal samples come from: ", paste(studies_normal, collapse = ", "))
  message("  This is a cross-cohort comparison. Interpret with caution.")
}

# ── Save main analysis data ──
main_expr <- expr[match(main_meta$sample, meta$sample), , drop = FALSE]

# Remove zero-variance genes
gene_vars <- apply(main_expr, 2, var, na.rm = TRUE)
zero_var <- which(gene_vars == 0)
if (length(zero_var) > 0) {
  message("Removing ", length(zero_var), " zero-variance genes")
  main_expr <- main_expr[, -zero_var, drop = FALSE]
}

saveRDS(main_meta, file.path(repo_root, "data", "processed", "metadata.rds"))
saveRDS(main_expr, file.path(repo_root, "data", "processed", "expression_matrix.rds"))

# ── Save TCGA-normal for sensitivity if available ──
if (n_tcga_norm > 0) {
  tcga_meta <- meta |> filter(proposed_condition %in% c("KIRP", "TCGA_Normal"))
  tcga_expr <- expr[match(tcga_meta$sample, meta$sample), , drop = FALSE]
  saveRDS(tcga_meta, file.path(repo_root, "data", "processed", "metadata_tcga_only.rds"))
  saveRDS(tcga_expr, file.path(repo_root, "data", "processed", "expression_matrix_tcga_only.rds"))
  message("TCGA-only sensitivity analysis data saved (", nrow(tcga_meta), " samples)")
}

message("\n✓ Sample preparation complete.")
message("  Main: ", nrow(main_meta), " samples (", n_kirp, " KIRP + ", n_gtx_norm, " Normal GTEx)")
if (n_tcga_norm > 0) message("  Sensitivity: TCGA only (", n_kirp, " KIRP + ", n_tcga_norm, " TCGA Normal)")
