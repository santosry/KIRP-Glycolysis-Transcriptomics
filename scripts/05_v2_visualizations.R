# 05_v2_visualizations.R
# v2.0.0 FIGURES AND TABLES — Central Carbon Metabolism

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(rio)
  library(ggplot2); library(ggrepel); library(pheatmap); library(RColorBrewer)
  library(tidyr)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
dir.create(file.path(repo_root, "results", "v2", "figures"), recursive = TRUE, showWarnings = FALSE)

# Load results
robustness <- rio::import(file.path(repo_root, "results", "v2", "robustness_classification.csv"))
deg_paired <- rio::import(file.path(repo_root, "results", "v2", "DEG_paired.csv"))
deg_gtx <- rio::import(file.path(repo_root, "results", "v2", "DEG_gtx.csv"))
deg_tcga <- rio::import(file.path(repo_root, "results", "v2", "DEG_tcga_all.csv"))
pw_map <- rio::import(file.path(repo_root, "results", "tables", "pathway_gene_membership.csv"))
meta <- rio::import(file.path(repo_root, "results", "v2", "sample_groups.csv"))

# Load expression data for heatmap
expr <- readRDS(file.path(repo_root, "data", "processed", "expression_matrix.rds"))
meta_expr <- readRDS(file.path(repo_root, "data", "processed", "metadata.rds"))

# Colors
col_robust <- "#1B5E20"
col_partial <- "#FF8F00"
col_sensitive <- "#1565C0"
col_discordant <- "#B71C1C"
col_inconclusive <- "grey70"
col_up <- "#1F5BFF"
col_down <- "#8A2BE2"

# Merge pathway info
robustness <- robustness |>
  left_join(pw_map |> select(gene, in_hsa00010, in_hsa00030, in_hsa00020), by = c("gene_id" = "gene"))

robustness$pathway_label <- ""
robustness$pathway_label[robustness$in_hsa00010] <- "Glycolysis/Gluconeogenesis"
robustness$pathway_label[robustness$in_hsa00030 & !robustness$in_hsa00010] <- "PPP"
robustness$pathway_label[robustness$in_hsa00020 & !robustness$in_hsa00010] <- "TCA"

# === FIGURE 1: Volcano-style robustness plot ===
robustness$mean_logFC <- rowMeans(
  cbind(robustness$logFC_paired, robustness$logFC_gtx, robustness$logFC_tcga), na.rm = TRUE
)
robustness$min_FDR <- pmin(robustness$FDR_paired, robustness$FDR_gtx, robustness$FDR_tcga, na.rm = TRUE)
robustness$min_FDR[is.na(robustness$min_FDR)] <- 1

rob_colors <- c(
  ROBUST = col_robust,
  PARTIALLY_ROBUST = col_partial,
  SENSITIVE_TO_COMPARATOR = col_sensitive,
  DIRECTION_DISCORDANT = col_discordant,
  INCONCLUSIVE = col_inconclusive
)

# Top genes to label
top_genes <- robustness |>
  filter(robustness_class %in% c("ROBUST", "PARTIALLY_ROBUST")) |>
  arrange(desc(abs(mean_logFC))) |>
  head(25) |>
  pull(gene_id)

robustness$label <- ifelse(robustness$gene_id %in% top_genes, robustness$gene_id, "")

p_robust <- ggplot(robustness, aes(mean_logFC, -log10(min_FDR), color = robustness_class)) +
  geom_point(aes(size = n_sig), alpha = 0.75) +
  scale_color_manual(values = rob_colors, name = "Robustness") +
  scale_size_continuous(range = c(1, 4), name = "N Comparators") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 30, box.padding = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
  labs(title = "Robustness of Transcriptomic Alterations in Central Carbon Metabolism",
       subtitle = paste0(sum(robustness$robustness_class == "ROBUST"), " robust | ",
                         sum(robustness$robustness_class == "PARTIALLY_ROBUST"), " partially robust | ",
                         sum(robustness$robustness_class == "DIRECTION_DISCORDANT"), " discordant"),
       x = "Mean log2 Fold Change (across 3 comparators)",
       y = "-log10(min FDR)") +
  theme_classic(base_size = 13)
ggsave(file.path(repo_root, "results", "v2", "figures", "Robustness_volcano.png"),
       p_robust, width = 10, height = 8, dpi = 300)

# === FIGURE 2: Comparator comparison — forest plot ===
robust_genes <- robustness |> filter(robustness_class %in% c("ROBUST", "PARTIALLY_ROBUST")) |>
  arrange(desc(abs(mean_logFC))) |> head(20)

forest_data <- robust_genes |>
  select(gene_id, logFC_paired, logFC_gtx, logFC_tcga) |>
  pivot_longer(-gene_id, names_to = "comparator", values_to = "logFC") |>
  mutate(
    comparator = case_when(
      comparator == "logFC_paired" ~ "Paired (KIRP, n=32 pairs)",
      comparator == "logFC_gtx" ~ "GTEx Normal (n=28)",
      comparator == "logFC_tcga" ~ "TCGA Normal (n=129)"
    ),
    gene_id = factor(gene_id, levels = rev(robust_genes$gene_id))
  )

p_forest <- ggplot(forest_data, aes(logFC, gene_id, color = comparator, shape = comparator)) +
  geom_point(size = 3, alpha = 0.85) +
  scale_color_manual(values = c("Paired (KIRP, n=32 pairs)" = col_robust,
                                 "GTEx Normal (n=28)" = col_sensitive,
                                 "TCGA Normal (n=129)" = col_partial)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  labs(title = "Robust and Partially Robust Genes — Effect Sizes by Comparator",
       x = "log2 Fold Change", y = "", color = "Comparator", shape = "Comparator") +
  theme_classic(base_size = 12)
ggsave(file.path(repo_root, "results", "v2", "figures", "Forest_robust_genes.png"),
       p_forest, width = 10, height = 8, dpi = 300)

# === FIGURE 3: Pathway-level heatmap ===
robust_and_partial <- robustness |>
  filter(robustness_class %in% c("ROBUST", "PARTIALLY_ROBUST")) |>
  arrange(robustness_class, desc(abs(mean_logFC)))

heat_genes <- robust_and_partial$gene_id
heat_genes <- intersect(heat_genes, colnames(expr))

# Use the main expression matrix (KIRP vs GTEx)
expr_heat <- expr[, heat_genes, drop = FALSE]
expr_scaled <- scale(expr_heat)
expr_scaled <- t(expr_scaled)

anno <- data.frame(
  Condition = meta_expr$condition,
  row.names = meta_expr$sample
)
anno_colors <- list(Condition = c(Normal_GTEx = col_up, KIRP = col_down))

# Add robustness annotation
rob_anno <- data.frame(
  Robustness = robust_and_partial$robustness_class[match(heat_genes, robust_and_partial$gene_id)],
  row.names = heat_genes
)
rob_cols <- c(ROBUST = col_robust, PARTIALLY_ROBUST = col_partial)

png(file.path(repo_root, "results", "v2", "figures", "Heatmap_robust_genes.png"),
    width = 12, height = 9, units = "in", res = 300)
pheatmap(expr_scaled,
         annotation_col = anno,
         annotation_colors = anno_colors,
         annotation_row = rob_anno,
         annotation_colors_row = list(Robustness = rob_cols),
         show_colnames = FALSE,
         main = paste0("Robust & Partially Robust Genes (n=", length(heat_genes), ") — Central Carbon Metabolism"),
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         clustering_method = "ward.D2",
         fontsize_row = 9)
dev.off()

# === FIGURE 4: Comparator agreement matrix ===
# Jaccard of DEG sets between comparators

get_degs <- function(df) {
  df$gene_id[df$regulation != "NS"]
}

degs_paired <- get_degs(deg_paired)
degs_gtx <- get_degs(deg_gtx)
degs_tcga <- get_degs(deg_tcga)

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

agreement_df <- tibble(
  Comparison = c("Paired vs GTEx", "Paired vs TCGA", "GTEx vs TCGA"),
  Jaccard = c(jaccard(degs_paired, degs_gtx),
              jaccard(degs_paired, degs_tcga),
              jaccard(degs_gtx, degs_tcga)),
  Shared = c(length(intersect(degs_paired, degs_gtx)),
             length(intersect(degs_paired, degs_tcga)),
             length(intersect(degs_gtx, degs_tcga))),
  Total_Union = c(length(union(degs_paired, degs_gtx)),
                  length(union(degs_paired, degs_tcga)),
                  length(union(degs_gtx, degs_tcga)))
)

rio::export(agreement_df, file.path(repo_root, "results", "v2", "comparator_agreement.csv"))

p_agree <- ggplot(agreement_df, aes(Comparison, Jaccard, fill = Comparison)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = paste0(Shared, "/", Total_Union, " genes")), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Paired vs GTEx" = "#D4A017", "Paired vs TCGA" = col_robust, "GTEx vs TCGA" = col_partial)) +
  labs(title = "DEG Overlap Between Comparators",
       subtitle = "Jaccard index of differentially expressed gene sets",
       y = "Jaccard Index") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(repo_root, "results", "v2", "figures", "Comparator_agreement.png"),
       p_agree, width = 8, height = 5, dpi = 300)

# === FIGURE 5: Pathway membership overview ===
pw_summary <- robustness |>
  filter(!is.na(mean_logFC)) |>
  mutate(
    pathway_group = case_when(
      in_hsa00010 & in_hsa00030 & in_hsa00020 ~ "All three",
      in_hsa00010 & in_hsa00030 ~ "Glycolysis + PPP",
      in_hsa00010 & in_hsa00020 ~ "Glycolysis + TCA",
      in_hsa00010 ~ "Glycolysis only",
      in_hsa00030 ~ "PPP only",
      in_hsa00020 ~ "TCA only",
      TRUE ~ "Unmapped"
    )
  )

pw_counts <- pw_summary |>
  group_by(pathway_group, robustness_class) |>
  summarise(n = n(), .groups = "drop") |>
  filter(pathway_group != "Unmapped")

p_pw <- ggplot(pw_counts, aes(pathway_group, n, fill = robustness_class)) +
  geom_col(alpha = 0.85) +
  scale_fill_manual(values = rob_colors) +
  labs(title = "Robustness by Pathway Membership",
       x = "Pathway Group", y = "Number of Genes", fill = "Robustness") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(repo_root, "results", "v2", "figures", "Pathway_robustness.png"),
       p_pw, width = 9, height = 6, dpi = 300)

# === Save robustness table ===
robustness_out <- robustness |>
  select(gene_id, mean_logFC, n_sig, robustness_class, 
         logFC_paired, FDR_paired, logFC_gtx, FDR_gtx, logFC_tcga, FDR_tcga,
         in_hsa00010, in_hsa00030, in_hsa00020) |>
  arrange(robustness_class, desc(abs(mean_logFC)))

rio::export(robustness_out, file.path(repo_root, "results", "v2", "robustness_table.csv"))

message("\n✓ v2 visualizations complete.")
message("  17 ROBUST | 6 PARTIALLY_ROBUST | 11 SENSITIVE | 7 DISCORDANT | 25 INCONCLUSIVE")
