suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tibble)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
meta_file <- file.path(repo_root, "data", "processed", "metadata.rds")
expr_file <- file.path(repo_root, "data", "processed", "expression_matrix.rds")

if (!file.exists(meta_file) || !file.exists(expr_file)) {
  stop("Run scripts/02_prepare_data.R before differential expression.")
}

dir.create(file.path(repo_root, "results", "differential_expression"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

meta <- readRDS(meta_file)
expr <- readRDS(expr_file)

E <- t(expr)
storage.mode(E) <- "numeric"

design <- model.matrix(~ 0 + condition, data = meta)
colnames(design) <- levels(meta$condition)
stopifnot(ncol(E) == nrow(design))

fit <- lmFit(E, design)
contrast_matrix <- makeContrasts(KIRP_vs_Normal = KIRP - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

deg <- topTable(
  fit2,
  coef = "KIRP_vs_Normal",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
) |>
  rownames_to_column(var = "gene_symbol")

deg$regulation <- "NS"
deg$regulation[deg$adj.P.Val < 0.05 & deg$logFC > 1] <- "Up"
deg$regulation[deg$adj.P.Val < 0.05 & deg$logFC < -1] <- "Down"

rio::export(deg, file.path(repo_root, "results", "differential_expression", "DEG_KIRP_vs_Normal.csv"))

summary_table <- as.data.frame(table(deg$regulation))
colnames(summary_table) <- c("regulation", "n")
rio::export(summary_table, file.path(repo_root, "results", "tables", "deg_summary.csv"))

message("Saved differential expression results.")
