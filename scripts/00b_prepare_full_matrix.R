# 00b_prepare_full_matrix.R
# Process the FULL expression matrix with all 109 pathway genes
# Replaces the old 66-gene preparation for 3-pathway analysis

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
raw_file <- file.path(repo_root, "data", "raw", "kidney_full.tsv")

if (!file.exists(raw_file)) stop("Run scripts/download_full_matrix.py first to get kidney_full.tsv")

message("Reading full matrix...")
raw <- rio::import(raw_file)
message(sprintf("Imported: %d rows x %d cols", nrow(raw), ncol(raw)))

# Gene columns (all except 'sample')
gene_cols <- setdiff(names(raw), "sample")
message(sprintf("Gene columns: %d", length(gene_cols)))

# Verify values are numeric
raw_num <- raw
for (gc in gene_cols) {
  raw_num[[gc]] <- as.numeric(as.character(raw_num[[gc]]))
}

# Check missing values
nas <- sum(is.na(raw_num[, gene_cols]))
message(sprintf("Missing values: %d", nas))

# Check scale — should be log2(norm_count+1)
vals <- unlist(raw_num[, gene_cols])
vals <- vals[is.finite(vals)]
message(sprintf("Expression range: [%.2f, %.2f], median: %.2f", min(vals), max(vals), median(vals)))

# Build expression matrix (samples × genes)
expr_matrix <- as.data.frame(raw_num[, gene_cols])
rownames(expr_matrix) <- raw$sample

# ── Build metadata from sample IDs ──
message("\nBuilding metadata from sample IDs...")
meta <- tibble(
  sample = raw$sample,
  condition = case_when(
    grepl("^GTEX-", sample) ~ "Normal_GTEx",
    grepl("-01", sample, fixed = TRUE) ~ "KIRP",
    grepl("-11", sample, fixed = TRUE) ~ "TCGA_Normal",
    TRUE ~ "Other"
  ),
  study = case_when(
    grepl("^GTEX-", sample) ~ "GTEX",
    TRUE ~ "TCGA"
  )
)

# Split TCGA_Normal into matched/unmatched
meta$participant <- sub("^(TCGA-..-....).*", "\\1", meta$sample)
kirp_participants <- unique(meta$participant[meta$condition == "KIRP"])
meta$paired_normal <- meta$condition == "TCGA_Normal" & meta$participant %in% kirp_participants

# Old condition-based grouping for backward compatibility  
meta$condition <- factor(meta$condition, levels = c("Normal_GTEx", "KIRP", "TCGA_Normal", "Other"))

message("Sample counts:")
print(table(meta$condition, meta$study))

message(sprintf("\nKIRP: %d | GTEx Normal: %d | TCGA Normal: %d (paired: %d)",
                sum(meta$condition == "KIRP"),
                sum(meta$condition == "Normal_GTEx"),
                sum(meta$condition == "TCGA_Normal"),
                sum(meta$paired_normal)))

# Confound check
des <- model.matrix(~ condition + study, data = meta[meta$condition %in% c("Normal_GTEx", "KIRP"), ])
message(sprintf("\nDesign matrix rank (condition+study): %d / %d cols", qr(des)$rank, ncol(des)))
if (qr(des)$rank < ncol(des)) {
  message("WARNING: condition and study are perfectly confounded!")
}

# ── Save ──
dir.create(file.path(repo_root, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
saveRDS(expr_matrix, file.path(repo_root, "data", "processed", "expression_matrix_full.rds"))
saveRDS(meta, file.path(repo_root, "data", "processed", "metadata_full.rds"))

# ── Also save pathway-specific matrices ──
pw <- rio::import(file.path(repo_root, "results", "tables", "pathway_gene_membership.csv"))

for (pw_id in c("hsa00010", "hsa00030", "hsa00020")) {
  pw_genes <- pw$gene[pw[[paste0("in_", pw_id)]]]
  pw_genes <- intersect(pw_genes, gene_cols)
  
  expr_pw <- expr_matrix[, pw_genes, drop = FALSE]
  saveRDS(expr_pw, file.path(repo_root, "data", "processed", paste0("expression_", pw_id, ".rds")))
  message(sprintf("  %s: %d genes saved", pw_id, length(pw_genes)))
}

# Write metadata for each pathway gene set
pw_summary <- tibble(
  pathway = c("hsa00010", "hsa00030", "hsa00020"),
  name = c("Glycolysis/Gluconeogenesis", "Pentose Phosphate Pathway", "Citrate Cycle (TCA)"),
  genes_in_kegg = c(
    sum(pw$in_hsa00010),
    sum(pw$in_hsa00030),
    sum(pw$in_hsa00020)
  ),
  genes_in_matrix = c(
    sum(pw$in_hsa00010 & pw$gene %in% gene_cols),
    sum(pw$in_hsa00030 & pw$gene %in% gene_cols),
    sum(pw$in_hsa00020 & pw$gene %in% gene_cols)
  )
)
rio::export(pw_summary, file.path(repo_root, "results", "tables", "pathway_coverage.csv"))

message("\nFull matrix preparation complete!")
message(sprintf("Total genes: %d | Samples: %d", length(gene_cols), nrow(meta)))
