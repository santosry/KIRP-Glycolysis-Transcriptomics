suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
meta_file <- file.path(repo_root, "data", "processed", "metadata.rds")
expr_file <- file.path(repo_root, "data", "processed", "expression_matrix.rds")
genes_file <- file.path(repo_root, "data", "metadata", "Hsa_genes.csv")

if (!file.exists(meta_file) || !file.exists(expr_file)) {
  stop("Run scripts/02_prepare_data.R before PCA/UMAP.")
}

dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

meta <- readRDS(meta_file)
expr <- readRDS(expr_file)

# transpose: samples x genes → genes x samples
E <- t(expr)
storage.mode(E) <- "numeric"

# remove constant genes (zero variance)
gene_vars <- apply(E, 1, var, na.rm = TRUE)
E_var <- E[gene_vars > 0, ]

message("Genes with variance > 0: ", nrow(E_var), " / ", nrow(E))

# ── PCA: all variable genes ──
pca_all <- prcomp(t(E_var), center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca_all$x)
pca_scores$condition <- meta$condition
pca_scores$sample <- meta$sample

var_explained <- round(100 * pca_all$sdev^2 / sum(pca_all$sdev^2), 1)

p_pca <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c(Normal = "#1F5BFF", KIRP = "#8A2BE2")) +
  labs(
    title = "PCA — Todos os genes com variância > 0",
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    color = "Grupo"
  ) +
  theme_classic(base_size = 14)

ggsave(
  filename = file.path(repo_root, "results", "figures", "PCA_all_genes.png"),
  plot = p_pca,
  width = 8,
  height = 6,
  dpi = 300
)

# ── PCA: only glycolysis/gluconeogenesis genes ──
if (file.exists(genes_file)) {
  pathway_genes <- rio::import(genes_file)$gene_symbol
  common_genes <- intersect(pathway_genes, rownames(E_var))

  if (length(common_genes) >= 3) {
    E_glyc <- E_var[common_genes, , drop = FALSE]
    pca_glyc <- prcomp(t(E_glyc), center = TRUE, scale. = TRUE)
    pca_glyc_scores <- as.data.frame(pca_glyc$x)
    pca_glyc_scores$condition <- meta$condition
    pca_glyc_scores$sample <- meta$sample

    var_glyc <- round(100 * pca_glyc$sdev^2 / sum(pca_glyc$sdev^2), 1)

    p_pca_glyc <- ggplot(pca_glyc_scores, aes(x = PC1, y = PC2, color = condition)) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_manual(values = c(Normal = "#1F5BFF", KIRP = "#8A2BE2")) +
      labs(
        title = paste0("PCA — Genes da via glicólise/gliconeogênese (n = ", length(common_genes), ")"),
        x = paste0("PC1 (", var_glyc[1], "%)"),
        y = paste0("PC2 (", var_glyc[2], "%)"),
        color = "Grupo"
      ) +
      theme_classic(base_size = 14)

    ggsave(
      filename = file.path(repo_root, "results", "figures", "PCA_glycolysis_genes.png"),
      plot = p_pca_glyc,
      width = 8,
      height = 6,
      dpi = 300
    )

    # PCA variance table
    pca_var_table <- tibble(
      PC = seq_along(var_glyc),
      variance_explained_pct = var_glyc
    )
    rio::export(pca_var_table, file.path(repo_root, "results", "tables", "pca_glycolysis_variance.csv"))

    message("Glycolysis PCA: ", length(common_genes), " genes used")
  } else {
    message("Fewer than 3 glycolysis genes with variance > 0; skipping focused PCA.")
  }
}

message("Saved PCA figures.")
