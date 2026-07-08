suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tibble)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
meta_file <- file.path(repo_root, "data", "processed", "metadata.rds")
expr_file <- file.path(repo_root, "data", "processed", "expression_matrix.rds")
genes_file <- file.path(repo_root, "data", "metadata", "Hsa_genes.csv")

if (!file.exists(meta_file) || !file.exists(expr_file)) {
  stop("Run scripts/02_prepare_data.R before differential expression.")
}
if (!file.exists(genes_file)) {
  stop("Run scripts/01_download_data.R before differential expression.")
}

dir.create(file.path(repo_root, "results", "differential_expression"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

meta <- readRDS(meta_file)
expr <- readRDS(expr_file)

E <- t(expr)
storage.mode(E) <- "numeric"

# ── Restrict to KEGG hsa00010 genes ──
pathway_genes <- rio::import(genes_file)$gene_symbol
common_genes <- intersect(pathway_genes, rownames(E))
if (length(common_genes) == 0) {
  stop("No KEGG hsa00010 genes found in the expression matrix.")
}
E <- E[common_genes, , drop = FALSE]
message("Genes analysed (intersection KEGG hsa00010 ∩ expression matrix): ", nrow(E))

# ── Model ──
# The input matrix is log2(norm_count + 1), already normalized and log-transformed
# as exported by UCSC Xena. In this scale, limma with eBayes() is appropriate:
# the log transformation stabilizes the mean-variance relationship, and the
# empirical Bayes moderation handles the small sample size in the normal group.
# For raw counts, voom() or limma-trend would be required (Ritchie et al., 2015).

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

n_up   <- sum(deg$regulation == "Up")
n_down <- sum(deg$regulation == "Down")
n_ns   <- sum(deg$regulation == "NS")
message(sprintf("DEGs: %d Up | %d Down | %d NS (total: %d genes analysed)", n_up, n_down, n_ns, nrow(deg)))
message("Saved differential expression results.")
