# 06_hsa00010_targeted_analysis.R
# EIXO 2: ANÁLISE DIRIGIDA DA VIA GLICÓLISE/GLICONEOGÊNESE

suppressPackageStartupMessages({
  library(KEGGREST)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(rio)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
deg_file <- file.path(repo_root, "results", "differential_expression", "DEG_global.csv")

if (!file.exists(deg_file)) {
  stop("Run 05_differential_expression_global.R first.")
}

dir.create(file.path(repo_root, "data", "metadata"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)

deg_all <- rio::import(deg_file)

# ── Download KEGG hsa00010 ──
message("Querying KEGG hsa00010...")
pw <- keggGet("hsa00010")
stopifnot(length(pw) == 1L)

genes_raw <- pw[[1]]$GENE
gene_symbols <- genes_raw[seq(2, length(genes_raw), by = 2)] |>
  str_remove("\\s*\\[.*$") |>
  str_trim() |>
  str_extract("^[^;]+") |>
  str_trim() |>
  unique() |>
  sort()

message("KEGG hsa00010 returned ", length(gene_symbols), " unique gene symbols")

# ── Save KEGG gene list ──
kegg_genes <- tibble(
  pathway_id = "hsa00010",
  gene_symbol = gene_symbols
)
rio::export(kegg_genes, file.path(repo_root, "data", "metadata", "Hsa_genes.csv"))

# ── Intersection with tested genes ──
tested_genes <- deg_all$gene_id
common <- intersect(gene_symbols, tested_genes)
not_found <- setdiff(gene_symbols, tested_genes)

message("Genes in hsa00010: ", length(gene_symbols))
message("Genes tested in model: ", length(tested_genes))
message("Intersection: ", length(common))
message("Not found in matrix: ", length(not_found))
if (length(not_found) > 0) message("  Missing: ", paste(not_found, collapse = ", "))

# ── Audit table for hsa00010 genes ──
audit <- tibble(
  gene_symbol = gene_symbols,
  in_matrix = gene_symbols %in% tested_genes,
  tested = gene_symbols %in% common
)

# Merge with DE results
audit <- audit |>
  left_join(deg_all |> select(gene_id, logFC, AveExpr, adj.P.Val, regulation),
            by = c("gene_symbol" = "gene_id"))

rio::export(audit, file.path(repo_root, "results", "tables", "hsa00010_audit.csv"))

# ── hsa00010 DEG counts ──
hsa_deg <- audit |> filter(tested & regulation %in% c("Up", "Down"))
n_up_hsa   <- sum(hsa_deg$regulation == "Up", na.rm = TRUE)
n_down_hsa <- sum(hsa_deg$regulation == "Down", na.rm = TRUE)
message("hsa00010 DEGs: ", n_up_hsa, " Up | ", n_down_hsa, " Down | ", nrow(hsa_deg), " total")

# ── Volcano plot highlighting hsa00010 ──
deg_all$in_hsa00010 <- deg_all$gene_id %in% common
deg_all$label <- ifelse(deg_all$gene_id %in% common & deg_all$regulation %in% c("Up", "Down"),
                         deg_all$gene_id, "")

deg_all$color <- deg_all$regulation
deg_all$color[deg_all$in_hsa00010 & deg_all$color == "NS"] <- "hsa_NS"
deg_all$color[deg_all$in_hsa00010 & deg_all$color == "Up"] <- "hsa_Up"
deg_all$color[deg_all$in_hsa00010 & deg_all$color == "Down"] <- "hsa_Down"

p_volcano_hsa <- ggplot(deg_all, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
  geom_point(size = 0.6, alpha = 0.5) +
  geom_point(data = filter(deg_all, in_hsa00010 & regulation %in% c("Up", "Down")),
             size = 2, alpha = 0.9, shape = 21, stroke = 0.8,
             fill = NA, color = "black") +
  scale_color_manual(values = c(
    Up = "#1F5BFF", Down = "#8A2BE2", NS = "grey80",
    hsa_Up = "#1F5BFF", hsa_Down = "#8A2BE2", hsa_NS = "grey60"
  )) +
  geom_text_repel(data = filter(deg_all, gene_id %in% common & regulation %in% c("Up", "Down")),
                  aes(label = gene_id), size = 3, max.overlaps = 50,
                  box.padding = 0.5, point.padding = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4) +
  labs(title = "Volcano plot global — Genes hsa00010 destacados",
       subtitle = paste0(n_up_hsa, " Up | ", n_down_hsa, " Down na via glicólise/gliconeogênese"),
       x = "log2 Fold Change", y = "-log10(FDR)") +
  guides(color = "none") +
  theme_classic(base_size = 14)

ggsave(file.path(repo_root, "results", "figures", "Volcano_hsa00010_highlighted.png"),
       p_volcano_hsa, width = 10, height = 8, dpi = 300)

# ── Heatmap of hsa00010 DEGs ──
if (nrow(hsa_deg) >= 3) {
  hsa_deg_genes <- hsa_deg$gene_symbol[!is.na(hsa_deg$gene_symbol)]
  
  expr_file <- file.path(repo_root, "data", "processed", "expression_matrix.rds")
  meta_file <- file.path(repo_root, "data", "processed", "metadata.rds")
  
  if (file.exists(expr_file) && file.exists(meta_file)) {
    expr_mat <- readRDS(expr_file)
    meta <- readRDS(meta_file)
    
    common_genes_expr <- intersect(hsa_deg_genes, colnames(expr_mat))
    if (length(common_genes_expr) >= 3) {
      expr_sub <- expr_mat[, common_genes_expr, drop = FALSE]
      
      # Scale genes (z-score)
      expr_scaled <- scale(expr_sub)
      expr_scaled <- t(expr_scaled)
      
      anno <- data.frame(Condition = meta$condition, row.names = meta$sample)
      anno_colors <- list(Condition = c(Normal = "#1F5BFF", KIRP = "#8A2BE2"))
      
      png(file.path(repo_root, "results", "figures", "Heatmap_hsa00010.png"),
          width = 12, height = 9, units = "in", res = 300)
      pheatmap(expr_scaled,
               annotation_col = anno,
               annotation_colors = anno_colors,
               show_colnames = FALSE,
               main = paste0("Genes hsa00010 DEGs (n = ", length(common_genes_expr), ")"),
               color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
               clustering_method = "ward.D2",
               fontsize_row = 8)
      dev.off()
      message("Heatmap saved: ", length(common_genes_expr), " genes")
    }
  }
}

# ── Ranking: hsa00010 genes by |logFC| ──
ranking <- audit |> filter(tested) |>
  mutate(abs_logFC = abs(logFC)) |>
  arrange(desc(abs_logFC))

rio::export(ranking, file.path(repo_root, "results", "tables", "hsa00010_ranking.csv"))

# Top 10
message("\nTop 10 hsa00010 genes by |logFC|:")
ranking_top <- head(ranking, 10)
for (i in 1:nrow(ranking_top)) {
  r <- ranking_top[i, ]
  message(sprintf("  %2d. %-10s logFC=% 6.2f  FDR=%.2e  %s",
                  i, r$gene_symbol, r$logFC, r$adj.P.Val, r$regulation))
}

message("\n✓ hsa00010 targeted analysis complete.")
