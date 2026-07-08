# 04_pca_umap.R
# PCA + UMAP com estabilidade quantitativa (Procrustes)

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(rio); library(ggplot2)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
meta_file <- file.path(repo_root, "data", "processed", "metadata.rds")
expr_file <- file.path(repo_root, "data", "processed", "expression_matrix.rds")
if (!file.exists(meta_file) || !file.exists(expr_file)) stop("Run 02_prepare_data.R first.")

dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

meta <- readRDS(meta_file); expr <- readRDS(expr_file)
E <- t(expr); storage.mode(E) <- "numeric"

# Pre-PCA missing check
if (anyNA(E)) stop("NA values in expression matrix. Cannot run PCA. Check data source.")

gene_vars <- apply(E, 1, var, na.rm = TRUE)
E_var <- E[gene_vars > 0, ]
n_genes_pca <- nrow(E_var)
message("Genes for PCA: ", n_genes_pca)

# ── PCA: scale=TRUE (genes standardized) ──
pca <- prcomp(t(E_var), center = TRUE, scale. = TRUE)
var_pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
message("PCA (scale=TRUE) | PC1: ", var_pct[1], "% | PC2: ", var_pct[2], "% | PC3: ", var_pct[3], "%")

# PCA scores
pca_scores <- as.data.frame(pca$x)
pca_scores$condition <- meta$condition
pca_scores$study <- meta$study
pca_scores$sample <- meta$sample
rio::export(pca_scores, file.path(repo_root, "results", "tables", "pca_scores.csv"))

# PCA loadings - top contributors
loadings <- as.data.frame(pca$rotation[, 1:min(5, ncol(pca$rotation))])
loadings$gene <- rownames(loadings)
rio::export(loadings, file.path(repo_root, "results", "tables", "pca_loadings.csv"))

top_pc1 <- head(loadings[order(-abs(loadings$PC1)), c("gene", "PC1")], 10)
message("Top PC1 genes: ", paste(top_pc1$gene, collapse = ", "))

# ── PCA plots (nano banana palette) ──
bn_bg    <- "#FFF8DC"; bn_panel <- "#FFFAF0"; bn_grid <- "#E8D5A3"
blue_n   <- "#1F5BFF"; purple_k <- "#8A2BE2"

p_pca <- ggplot(pca_scores, aes(PC1, PC2, color = condition)) +
  geom_point(size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c(Normal_GTEx = blue_n, KIRP = purple_k)) +
  labs(title = "PCA — Global Transcriptome",
       subtitle = paste0("PC1 (", var_pct[1], "%) vs PC2 (", var_pct[2], "%) | ", n_genes_pca, " genes"),
       x = paste0("PC1 (", var_pct[1], "%)"), y = paste0("PC2 (", var_pct[2], "%)")) +
  theme_minimal(base_size = 14) +
  theme(plot.background = element_rect(fill = bn_bg, color = NA),
        panel.background = element_rect(fill = bn_panel),
        panel.grid = element_line(color = bn_grid, linewidth = 0.3))
ggsave(file.path(repo_root, "results", "figures", "PCA_PC1_PC2.png"), p_pca, width = 8, height = 6, dpi = 300)

# PCA by study
p_pca_s <- ggplot(pca_scores, aes(PC1, PC2, color = study)) +
  geom_point(size = 2.5, alpha = 0.75) + scale_color_brewer(palette = "Set1") +
  labs(title = "PCA — Colored by Study Cohort", x = paste0("PC1 (", var_pct[1], "%)"),
       y = paste0("PC2 (", var_pct[2], "%)")) +
  theme_minimal(base_size = 14) +
  theme(plot.background = element_rect(fill = bn_bg, color = NA),
        panel.background = element_rect(fill = bn_panel),
        panel.grid = element_line(color = bn_grid, linewidth = 0.3))
ggsave(file.path(repo_root, "results", "figures", "PCA_by_study.png"), p_pca_s, width = 8, height = 6, dpi = 300)

# Scree
scree_df <- tibble(PC = seq_along(var_pct), Variance = var_pct, Cumulative = cumsum(var_pct))
n_show <- min(20, nrow(scree_df))
p_scree <- ggplot(scree_df[1:n_show, ], aes(PC, Variance)) +
  geom_col(fill = "#D4A017", alpha = 0.8) +
  geom_line(aes(y = Cumulative), color = purple_k, linewidth = 1) +
  geom_point(aes(y = Cumulative), color = purple_k, size = 2) +
  labs(title = "Scree Plot", x = "Principal Component", y = "Variance Explained (%)") +
  scale_x_continuous(breaks = 1:n_show) + theme_minimal(base_size = 13)
ggsave(file.path(repo_root, "results", "figures", "PCA_scree.png"), p_scree, width = 10, height = 5, dpi = 200)
rio::export(scree_df, file.path(repo_root, "results", "tables", "pca_variance.csv"))

# ── UMAP with quantitative stability ──
if (!requireNamespace("umap", quietly = TRUE)) {
  message("umap package not available — skipping UMAP")
} else {
  library(umap)
  
  # Use top PCs as input (more stable than raw genes)
  n_pcs_umap <- min(50, ncol(pca$x))
  umap_input <- pca$x[, 1:n_pcs_umap]
  message("UMAP input: ", n_pcs_umap, " PCs")
  
  # Stability: same params, multiple seeds
  set.seed(1); u1 <- umap(umap_input, n_neighbors = 15, min_dist = 0.1)$layout
  set.seed(42); u2 <- umap(umap_input, n_neighbors = 15, min_dist = 0.1)$layout
  set.seed(123); u3 <- umap(umap_input, n_neighbors = 15, min_dist = 0.1)$layout
  
  # Procrustes analysis between embeddings
  procrustes_corr <- function(X, Y) {
    # Align Y to X via Procrustes, return correlation
    if (!requireNamespace("vegan", quietly = TRUE)) {
      # Simple: just correlate distance matrices
      dx <- as.matrix(dist(X)); dy <- as.matrix(dist(Y))
      return(cor(dx[lower.tri(dx)], dy[lower.tri(dy)]))
    }
    proc <- vegan::procrustes(X, Y, scale = FALSE)
    return(cor(c(X), c(fitted(proc))))
  }
  
  r12 <- procrustes_corr(u1, u2)
  r13 <- procrustes_corr(u1, u3)
  r23 <- procrustes_corr(u2, u3)
  message(sprintf("UMAP Procrustes correlations: r12=%.3f r13=%.3f r23=%.3f", r12, r13, r23))
  
  stab_df <- tibble(
    metric = c("n_pcs_input", "n_neighbors", "min_dist", "procrustes_r12", "procrustes_r13", "procrustes_r23"),
    value = c(n_pcs_umap, 15, 0.1, round(r12, 4), round(r13, 4), round(r23, 4))
  )
  rio::export(stab_df, file.path(repo_root, "results", "tables", "umap_stability.csv"))
  
  # Main UMAP (seed 1)
  umap_df <- as.data.frame(u1); colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$condition <- meta$condition
  umap_df$study <- meta$study
  umap_df$sample <- meta$sample
  rio::export(umap_df, file.path(repo_root, "results", "tables", "umap_coordinates.csv"))
  
  # UMAP plots
  p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2, color = condition)) +
    geom_point(size = 2.5, alpha = 0.75) +
    scale_color_manual(values = c(Normal_GTEx = blue_n, KIRP = purple_k)) +
    labs(title = "UMAP — Global Transcriptome", subtitle = paste0(n_pcs_umap, " PCs | n_neighbors=15 min_dist=0.1 | Procrustes r=", round(mean(c(r12,r13,r23)), 3)),
         color = "Group") +
    theme_minimal(base_size = 14) +
    theme(plot.background = element_rect(fill = bn_bg, color = NA),
          panel.background = element_rect(fill = bn_panel),
          panel.grid = element_line(color = bn_grid, linewidth = 0.3))
  ggsave(file.path(repo_root, "results", "figures", "UMAP.png"), p_umap, width = 8, height = 6, dpi = 300)
  
  p_umap_s <- ggplot(umap_df, aes(UMAP1, UMAP2, color = study)) +
    geom_point(size = 2.5, alpha = 0.75) + scale_color_brewer(palette = "Set1") +
    labs(title = "UMAP — Colored by Study Cohort") +
    theme_minimal(base_size = 14) +
    theme(plot.background = element_rect(fill = bn_bg, color = NA),
          panel.background = element_rect(fill = bn_panel),
          panel.grid = element_line(color = bn_grid, linewidth = 0.3))
  ggsave(file.path(repo_root, "results", "figures", "UMAP_by_study.png"), p_umap_s, width = 8, height = 6, dpi = 300)
  
  # Confounding note
  studies_k <- unique(meta$study[meta$condition == "KIRP"])
  studies_n <- unique(meta$study[meta$condition == "Normal_GTEx"])
  if (length(intersect(studies_k, studies_n)) == 0) {
    message("NOTE: condition and study perfectly confounded. UMAP colorings by condition and study represent the same structural separation.")
  }
}

message("\n✓ PCA & UMAP complete. Genes: ", n_genes_pca, " | PC1: ", var_pct[1], "% | PC2: ", var_pct[2], "%")
