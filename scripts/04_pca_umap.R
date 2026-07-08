# 04_pca_umap.R
# PCA + UMAP com análise de estabilidade e confundimento TCGA vs GTEx

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
meta_file <- file.path(repo_root, "data", "processed", "metadata.rds")
expr_file <- file.path(repo_root, "data", "processed", "expression_matrix.rds")

if (!file.exists(meta_file) || !file.exists(expr_file)) {
  stop("Run 02_prepare_data.R first.")
}

dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

meta <- readRDS(meta_file)
expr <- readRDS(expr_file)

# Transpose: samples × genes → genes × samples
E <- t(expr)
storage.mode(E) <- "numeric"

# Remove constant genes
gene_vars <- apply(E, 1, var, na.rm = TRUE)
E_var <- E[gene_vars > 0, ]
n_genes_used <- nrow(E_var)
message("Genes with variance > 0: ", n_genes_used)

# ── PCA ──
pca <- prcomp(t(E_var), center = TRUE, scale. = TRUE)
var_pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

message("Variance explained:")
message("  PC1: ", var_pct[1], "%")
message("  PC2: ", var_pct[2], "%")
message("  PC3: ", var_pct[3], "%")
message("  Cumulative PC1-3: ", sum(var_pct[1:3]), "%")

pca_scores <- as.data.frame(pca$x)
pca_scores$condition <- meta$condition
pca_scores$study <- meta$study
pca_scores$sample <- meta$sample

# ── 🍌 Nano Banana Palette ──
yellow_light  <- "#FFF8DC"
yellow_mid    <- "#FFE135"
yellow_dark   <- "#D4A017"
purple_kirp   <- "#8A2BE2"
blue_normal   <- "#1F5BFF"
grey_bg       <- "#F5F0E1"

# ── PCA PC1 × PC2 ──
p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 2.5, alpha = 0.75) +
  scale_color_manual(values = c(Normal = blue_normal, KIRP = purple_kirp)) +
  labs(
    title = "PCA — Transcriptoma global",
    subtitle = paste0("PC1 (", var_pct[1], "%) vs PC2 (", var_pct[2], "%) | ", n_genes_used, " genes"),
    x = paste0("PC1 (", var_pct[1], "%)"),
    y = paste0("PC2 (", var_pct[2], "%)"),
    color = "Grupo"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = grey_bg, color = NA),
    panel.background = element_rect(fill = "#FFFAF0"),
    panel.grid = element_line(color = "#E8D5A3", linewidth = 0.3)
  )

ggsave(file.path(repo_root, "results", "figures", "PCA_PC1_PC2.png"), p_pca, width = 8, height = 6, dpi = 300)

# ── PCA colored by study ──
p_pca_study <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = study)) +
  geom_point(size = 2.5, alpha = 0.75) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "PCA — Colorido por coorte de origem (study)",
       subtitle = "Verificar confundimento TCGA × GTEx",
       x = paste0("PC1 (", var_pct[1], "%)"),
       y = paste0("PC2 (", var_pct[2], "%)"),
       color = "Study") +
  theme_minimal(base_size = 14) +
  theme(
    plot.background = element_rect(fill = grey_bg, color = NA),
    panel.background = element_rect(fill = "#FFFAF0"),
    panel.grid = element_line(color = "#E8D5A3", linewidth = 0.3)
  )

ggsave(file.path(repo_root, "results", "figures", "PCA_by_study.png"), p_pca_study, width = 8, height = 6, dpi = 300)

# ── Scree plot ──
scree_data <- tibble(PC = seq_along(var_pct), Variance = var_pct)
n_show <- min(20, length(var_pct))

p_scree <- ggplot(scree_data[1:n_show, ], aes(x = PC, y = Variance)) +
  geom_col(fill = yellow_dark, alpha = 0.8) +
  geom_line(aes(y = cumsum(Variance)), color = purple_kirp, linewidth = 1) +
  geom_point(aes(y = cumsum(Variance)), color = purple_kirp, size = 2) +
  labs(title = "Scree plot — Variância explicada por PC",
       x = "Componente Principal", y = "Variância explicada (%)") +
  scale_x_continuous(breaks = 1:n_show) +
  theme_minimal(base_size = 13) +
  theme(plot.background = element_rect(fill = grey_bg, color = NA))

ggsave(file.path(repo_root, "results", "figures", "PCA_scree.png"), p_scree, width = 10, height = 5, dpi = 200)

# ── PCA variance table ──
pca_var_table <- tibble(PC = seq_along(var_pct), variance_pct = var_pct,
                         cumulative_pct = cumsum(var_pct))
rio::export(pca_var_table, file.path(repo_root, "results", "tables", "pca_variance.csv"))

# ── UMAP ──
if (requireNamespace("umap", quietly = TRUE)) {
  umap_configs <- list(
    list(seed = 1,  n_neighbors = 15, min_dist = 0.1),
    list(seed = 42, n_neighbors = 15, min_dist = 0.1),
    list(seed = 1,  n_neighbors = 30, min_dist = 0.3)
  )
  
  for (cfg in umap_configs) {
    set.seed(cfg$seed)
    umap_out <- umap::umap(t(E_var),
                           n_neighbors = cfg$n_neighbors,
                           min_dist = cfg$min_dist,
                           metric = "euclidean")
    
    umap_df <- as.data.frame(umap_out$layout)
    colnames(umap_df) <- c("UMAP1", "UMAP2")
    umap_df$condition <- meta$condition
    umap_df$study <- meta$study
    umap_df$sample <- meta$sample
    
    p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = condition)) +
      geom_point(size = 2.5, alpha = 0.75) +
      scale_color_manual(values = c(Normal = blue_normal, KIRP = purple_kirp)) +
      labs(title = paste0("UMAP — Transcriptoma global"),
           subtitle = paste0("seed=", cfg$seed, " | n_neighbors=", cfg$n_neighbors,
                           " | min_dist=", cfg$min_dist, " | ", n_genes_used, " genes"),
           color = "Grupo") +
      theme_minimal(base_size = 14) +
      theme(
        plot.background = element_rect(fill = grey_bg, color = NA),
        panel.background = element_rect(fill = "#FFFAF0"),
        panel.grid = element_line(color = "#E8D5A3", linewidth = 0.3)
      )
    
    suffix <- paste0("s", cfg$seed, "_n", cfg$n_neighbors, "_d", cfg$min_dist*10)
    ggsave(file.path(repo_root, "results", "figures", paste0("UMAP_", suffix, ".png")),
           p_umap, width = 8, height = 6, dpi = 300)
    
    message("UMAP saved: ", suffix)
  }
} else {
  message("umap package not available — skipping UMAP")
}

# ── Confounding audit ──
confound <- as.data.frame(table(meta$condition, meta$study))
colnames(confound) <- c("condition", "study", "n")
rio::export(confound, file.path(repo_root, "results", "tables", "confounding_audit.csv"))

message("\n✓ PCA & UMAP complete. Genes used: ", n_genes_used)
message("  PC1: ", var_pct[1], "% | PC2: ", var_pct[2], "%")
