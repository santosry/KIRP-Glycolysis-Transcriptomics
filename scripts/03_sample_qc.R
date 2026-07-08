# 03_sample_qc.R
# QC DAS AMOSTRAS: boxplots, density, correlation, clustering

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
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

# expression matrix: samples × genes
message("Samples: ", nrow(expr), " × Genes: ", ncol(expr))

# ── Per-sample QC metrics ──
qc <- tibble(
  sample     = meta$sample,
  condition  = meta$condition,
  study      = meta$study,
  median_expr = apply(expr, 1, median, na.rm = TRUE),
  iqr_expr    = apply(expr, 1, IQR, na.rm = TRUE),
  var_expr    = apply(expr, 1, var, na.rm = TRUE),
  min_expr    = apply(expr, 1, min, na.rm = TRUE),
  max_expr    = apply(expr, 1, max, na.rm = TRUE),
  pct_zero    = apply(expr, 1, function(x) mean(x == 0, na.rm = TRUE) * 100)
)

# Flag potential outliers: samples with median > 2 SD from group median
qc <- qc |> group_by(condition) |>
  mutate(
    median_z = (median_expr - mean(median_expr)) / sd(median_expr),
    flag_outlier = abs(median_z) > 3
  ) |> ungroup()

n_outliers <- sum(qc$flag_outlier)
if (n_outliers > 0) {
  message("Potential outliers flagged (>3 SD): ", n_outliers)
  message("  Samples: ", paste(qc$sample[qc$flag_outlier], collapse = ", "))
}

rio::export(qc, file.path(repo_root, "results", "tables", "sample_qc.csv"))

# ── Boxplot: expression distribution by sample ──
# downsample for plotting if too many samples
sample_order <- order(qc$median_expr)
expr_melt <- reshape2::melt(expr)
colnames(expr_melt) <- c("sample_idx", "gene_idx", "expression")
expr_melt$sample <- meta$sample[expr_melt$sample_idx]
expr_melt$condition <- meta$condition[expr_melt$sample_idx]

p_box <- ggplot(expr_melt, aes(x = reorder(sample, expression, median), y = expression, fill = condition)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.1) +
  scale_fill_manual(values = c(Normal = "#1F5BFF", KIRP = "#8A2BE2")) +
  labs(x = "Amostra", y = "Expressão (log2)", fill = "Grupo",
       title = "Distribuição de expressão por amostra") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(file.path(repo_root, "results", "figures", "QC_boxplot.png"), p_box, width = 14, height = 5, dpi = 200)

# ── Density plot ──
p_dens <- ggplot(expr_melt, aes(x = expression, color = condition, group = sample)) +
  geom_density(alpha = 0.05, linewidth = 0.15) +
  scale_color_manual(values = c(Normal = "#1F5BFF", KIRP = "#8A2BE2")) +
  labs(x = "Expressão (log2)", y = "Densidade", color = "Grupo",
       title = "Distribuição de densidade por amostra") +
  theme_classic(base_size = 14)

ggsave(file.path(repo_root, "results", "figures", "QC_density.png"), p_dens, width = 8, height = 6, dpi = 200)

# ── Correlation heatmap ──
# Use top variable genes for correlation
top_genes <- order(apply(expr, 2, var, na.rm = TRUE), decreasing = TRUE)[1:min(1000, ncol(expr))]
cor_mat <- cor(t(expr[, top_genes]), use = "pairwise.complete.obs")

# Annotation
anno <- data.frame(
  Condition = meta$condition,
  Study = meta$study,
  row.names = meta$sample
)
anno_colors <- list(
  Condition = c(Normal = "#1F5BFF", KIRP = "#8A2BE2"),
  Study = c(TCGA = "#FF6B6B", GTEX = "#4ECDC4")
)

png(file.path(repo_root, "results", "figures", "QC_correlation_heatmap.png"),
    width = 10, height = 8, units = "in", res = 200)
pheatmap(cor_mat,
         annotation_col = anno,
         annotation_colors = anno_colors,
         show_rownames = FALSE, show_colnames = FALSE,
         main = "Correlação amostra-amostra (top 1000 genes variáveis)",
         clustering_method = "ward.D2",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
dev.off()

# ── QC summary table ──
qc_summary <- tibble(
  metric = c("n_samples_total", "n_KIRP", "n_Normal",
             "n_outliers_flagged", "pct_outliers",
             "median_expression_KIRP", "median_expression_Normal",
             "iqr_expression_KIRP", "iqr_expression_Normal"),
  value = c(
    as.character(nrow(meta)), as.character(sum(meta$condition == "KIRP")),
    as.character(sum(meta$condition == "Normal")),
    as.character(n_outliers), as.character(round(100*n_outliers/nrow(meta), 1)),
    as.character(round(median(qc$median_expr[qc$condition=="KIRP"]), 2)),
    as.character(round(median(qc$median_expr[qc$condition=="Normal"]), 2)),
    as.character(round(median(qc$iqr_expr[qc$condition=="KIRP"]), 2)),
    as.character(round(median(qc$iqr_expr[qc$condition=="Normal"]), 2))
  )
)
rio::export(qc_summary, file.path(repo_root, "results", "tables", "qc_summary.csv"))

message("\n✓ Sample QC complete.")
message("  KIRP samples: ", sum(meta$condition == "KIRP"))
message("  Normal samples: ", sum(meta$condition == "Normal"))
message("  Flagged outliers: ", n_outliers)
