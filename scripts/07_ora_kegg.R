# 07_ora_kegg.R
# ORA KEGG — separadamente para Up e Down (usando universo transcriptômico global)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(dplyr)
  library(tibble)
  library(rio)
  library(ggplot2)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
deg_file <- file.path(repo_root, "results", "differential_expression", "DEG_global.csv")
universe_file <- file.path(repo_root, "results", "differential_expression", "gene_universe.csv")

if (!file.exists(deg_file)) stop("Run 05_differential_expression_global.R first.")

dir.create(file.path(repo_root, "results", "enrichment"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

deg_all <- rio::import(deg_file)

# ── Gene ID conversion ──
# Input: HGNC symbols → Entrez
universe_genes <- readLines(universe_file)
universe_genes <- universe_genes[-1]  # remove header if present

# Try simple approach: use Bitr
message("Converting gene symbols to Entrez IDs...")

# Get universe Entrez
universe_entrez <- tryCatch({
  bitr(universe_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}, error = function(e) {
  message("bitr failed, trying alternative...")
  NULL
})

if (is.null(universe_entrez) || nrow(universe_entrez) == 0) {
  # Try with keys
  all_keys <- keys(org.Hs.eg.db, keytype = "SYMBOL")
  valid_genes <- intersect(universe_genes, all_keys)
  universe_entrez <- select(org.Hs.eg.db, keys = valid_genes,
                            columns = "ENTREZID", keytype = "SYMBOL")
  universe_entrez <- universe_entrez[!is.na(universe_entrez$ENTREZID), ]
  universe_entrez <- universe_entrez[!duplicated(universe_entrez$ENTREZID), ]
}

message("Universe Entrez IDs: ", nrow(universe_entrez))

# ── ORA KEGG: UP ──
up_genes <- deg_all$gene_id[deg_all$regulation == "Up"]
up_entrez <- universe_entrez$ENTREZID[universe_entrez$SYMBOL %in% up_genes]
message("Up DEGs with Entrez: ", length(up_entrez), " / ", length(up_genes))

kegg_up <- NULL
if (length(up_entrez) >= 5) {
  kegg_up <- enrichKEGG(
    gene = up_entrez,
    universe = universe_entrez$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
    kegg_up <- as_tibble(kegg_up@result)
    rio::export(kegg_up, file.path(repo_root, "results", "enrichment", "KEGG_ORA_Up.csv"))
    message("KEGG Up: ", nrow(kegg_up), " terms")
  }
} else {
  message("Too few Up DEGs for ORA (", length(up_entrez), ")")
}

# ── ORA KEGG: DOWN ──
down_genes <- deg_all$gene_id[deg_all$regulation == "Down"]
down_entrez <- universe_entrez$ENTREZID[universe_entrez$SYMBOL %in% down_genes]
message("Down DEGs with Entrez: ", length(down_entrez), " / ", length(down_genes))

kegg_down <- NULL
if (length(down_entrez) >= 5) {
  kegg_down <- enrichKEGG(
    gene = down_entrez,
    universe = universe_entrez$ENTREZID,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
    kegg_down <- as_tibble(kegg_down@result)
    rio::export(kegg_down, file.path(repo_root, "results", "enrichment", "KEGG_ORA_Down.csv"))
    message("KEGG Down: ", nrow(kegg_down), " terms")
  }
} else {
  message("Too few Down DEGs for ORA (", length(down_entrez), ")")
}

# ── Dotplots ──
plot_ora <- function(ora_df, direction, max_terms = 20) {
  if (is.null(ora_df) || nrow(ora_df) == 0) return(NULL)
  
  df <- ora_df |> filter(p.adjust < 0.05)
  if (nrow(df) == 0) {
    message("No significant terms for ", direction)
    return(NULL)
  }
  
  df <- head(df[order(df$p.adjust), ], max_terms)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  
  ggplot(df, aes(x = Count, y = Description, size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#8A2BE2", high = "#FFE135",
                         trans = "log10", name = "FDR") +
    labs(title = paste0("KEGG ORA — ", direction),
         x = "Gene Count", y = "") +
    theme_minimal(base_size = 12) +
    theme(plot.background = element_rect(fill = "#FFF8DC", color = NA))
}

if (!is.null(kegg_up)) {
  p_up <- plot_ora(kegg_up, "Upregulated")
  if (!is.null(p_up)) ggsave(file.path(repo_root, "results", "figures", "KEGG_dotplot_Up.png"),
                              p_up, width = 10, height = 7, dpi = 300)
}

if (!is.null(kegg_down)) {
  p_down <- plot_ora(kegg_down, "Downregulated")
  if (!is.null(p_down)) ggsave(file.path(repo_root, "results", "figures", "KEGG_dotplot_Down.png"),
                                p_down, width = 10, height = 7, dpi = 300)
}

# ── Check if glycolysis appears ──
glycolysis_ids <- c("hsa00010")
for (gid in glycolysis_ids) {
  in_up <- !is.null(kegg_up) && gid %in% kegg_up$ID
  in_down <- !is.null(kegg_down) && gid %in% kegg_down$ID
  message(sprintf("  %s | In Up: %s | In Down: %s | p.adjust Up: %s | p.adjust Down: %s",
                  gid, in_up, in_down,
                  if(in_up) format(kegg_up$p.adjust[kegg_up$ID == gid], digits=3) else "N/A",
                  if(in_down) format(kegg_down$p.adjust[kegg_down$ID == gid], digits=3) else "N/A"))
}

message("\n✓ KEGG ORA complete.")
