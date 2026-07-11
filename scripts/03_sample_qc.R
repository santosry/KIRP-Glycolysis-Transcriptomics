# 03_sample_qc.R
# QC DAS AMOSTRAS — sem reshape2 (usa tidyr::pivot_longer)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
meta_file <- file.path(repo_root, "data", "processed", "metadata.rds")
expr_file <- file.path(repo_root, "data", "processed", "expression_matrix.rds")
if (!file.exists(meta_file) || !file.exists(expr_file)) stop("Run 02_prepare_data.R first.")

dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

meta <- readRDS(meta_file)
expr <- readRDS(expr_file)
message("Samples: ", nrow(expr), " × Genes: ", ncol(expr))

# ── Per-sample QC ──
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

qc <- qc |> group_by(condition) |>
  mutate(median_z = (median_expr - mean(median_expr)) / sd(median_expr),
         flag_outlier = abs(median_z) > 3) |> ungroup()

n_outliers <- sum(qc$flag_outlier)
if (n_outliers > 0) message("Potential outliers (>3 SD): ", n_outliers)

rio::export(qc, file.path(repo_root, "results", "tables", "sample_qc.csv"))

# ── Boxplot — usando pivot_longer do tidyr ──
expr_long <- as.data.frame(expr) |>
  mutate(sample = meta$sample, condition = meta$condition) |>
  pivot_longer(cols = -c(sample, condition), names_to = "gene", values_to = "expression")

p_box <- ggplot(expr_long, aes(x = reorder(sample, expression, median), y = expression, fill = condition)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.1) +
  scale_fill_manual(values = c(Normal_GTEx = "#1F5BFF", KIRP = "#8A2BE2")) +
  labs(x = "Sample", y = "Expression (log2)", fill = "Group",
       title = "Expression distribution by sample") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(repo_root, "results", "figures", "QC_boxplot.png"), p_box, width = 14, height = 5, dpi = 200)

# ── Density ──
p_dens <- ggplot(expr_long, aes(x = expression, color = condition, group = sample)) +
  geom_density(alpha = 0.05, linewidth = 0.15) +
  scale_color_manual(values = c(Normal_GTEx = "#1F5BFF", KIRP = "#8A2BE2")) +
  labs(x = "Expression (log2)", y = "Density", color = "Group",
       title = "Density distribution by sample") +
  theme_classic(base_size = 14)
ggsave(file.path(repo_root, "results", "figures", "QC_density.png"), p_dens, width = 8, height = 6, dpi = 200)

# ── Correlation heatmap ──
top_genes <- order(apply(expr, 2, var, na.rm = TRUE), decreasing = TRUE)[1:min(1000, ncol(expr))]
cor_mat <- cor(t(expr[, top_genes]), use = "pairwise.complete.obs")

anno <- data.frame(Condition = meta$condition, Study = meta$study, row.names = meta$sample)
anno_colors <- list(Condition = c(Normal_GTEx = "#1F5BFF", KIRP = "#8A2BE2"),
                    Study = c(TCGA = "#FF6B6B", GTEX = "#4ECDC4"))

png(file.path(repo_root, "results", "figures", "QC_correlation_heatmap.png"), width = 10, height = 8, units = "in", res = 200)
pheatmap(cor_mat, annotation_col = anno, annotation_colors = anno_colors,
         show_rownames = FALSE, show_colnames = FALSE,
         main = "Sample correlation (top 1000 variable genes)",
         clustering_method = "ward.D2",
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
dev.off()

# ── Missing values check ──
n_missing <- sum(is.na(expr))
message("Missing values in expression matrix: ", n_missing)
if (n_missing > 0) {
  missing_by_gene <- colSums(is.na(expr))
  missing_by_sample <- rowSums(is.na(expr))
  rio::export(data.frame(gene = names(missing_by_gene), n_missing = missing_by_gene),
              file.path(repo_root, "results", "tables", "missing_by_gene.csv"))
  rio::export(data.frame(sample = meta$sample, n_missing = missing_by_sample),
              file.path(repo_root, "results", "tables", "missing_by_sample.csv"))
  stop("Missing values detected. Pipeline cannot proceed with NA in expression matrix. Check data source.")
}

# ── QC summary ──
qc_summary <- tibble(
  metric = c("n_samples_total", "n_KIRP", "n_Normal_GTEx",
             "n_outliers_flagged", "pct_outliers",
             "median_expression_KIRP", "median_expression_Normal",
             "n_missing_values"),
  value = c(nrow(meta), sum(meta$condition == "KIRP"), sum(meta$condition == "Normal_GTEx"),
            n_outliers, round(100*n_outliers/nrow(meta), 1),
            round(median(qc$median_expr[qc$condition=="KIRP"]), 2),
            round(median(qc$median_expr[qc$condition=="Normal_GTEx"]), 2),
            n_missing)
)
rio::export(qc_summary, file.path(repo_root, "results", "tables", "qc_summary.csv"))
message("\n✓ Sample QC complete.")
