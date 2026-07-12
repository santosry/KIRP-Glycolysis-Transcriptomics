# 09_rank_based_enrichment.R
# GSEA KEGG + Reactome usando moderated t-statistic

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA)
  library(enrichplot); library(dplyr); library(tibble); library(rio); library(ggplot2)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
deg_file <- file.path(repo_root, "results", "differential_expression", "DEG_global.csv")
mapping_file <- file.path(repo_root, "results", "tables", "gene_id_mapping.csv")

if (!file.exists(deg_file)) stop("Run 05_differential_expression_global.R first.")
if (!file.exists(mapping_file)) stop("Run 07_ora_kegg.R first (generates gene_id_mapping.csv).")

dir.create(file.path(repo_root, "results", "enrichment"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)

deg_all <- rio::import(deg_file)
mapping <- rio::import(mapping_file)

# Build ranked list: moderated t-statistic, ENTREZID
# Map gene_id → ENTREZID using the mapping
# mapping contains gene_symbol + ENTREZID (already cleaned)
deg_mapped <- deg_all |>
  inner_join(mapping, by = c("gene_id" = "gene_symbol"))

# Use moderated t (t column from limma) as ranking statistic
ranked_list <- deg_mapped$t
names(ranked_list) <- deg_mapped$ENTREZID
ranked_list <- sort(ranked_list, decreasing = TRUE)

message("GSEA input: ", length(ranked_list), " genes with moderated t-statistic")

# ── GSEA KEGG ──
gsea_kegg <- gseKEGG(
  geneList = ranked_list,
  organism = "hsa",
  pvalueCutoff = 1,
  eps = 0
)

if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  gsea_kegg_df <- as_tibble(gsea_kegg@result)
  rio::export(gsea_kegg_df, file.path(repo_root, "results", "enrichment", "GSEA_KEGG.csv"))
  message("GSEA KEGG: ", nrow(gsea_kegg_df), " gene sets")
  
  # Check hsa00010
  if ("hsa00010" %in% gsea_kegg_df$ID) {
    message("  hsa00010 NES: ", round(gsea_kegg_df$NES[gsea_kegg_df$ID=="hsa00010"], 3),
            " | FDR: ", format(gsea_kegg_df$p.adjust[gsea_kegg_df$ID=="hsa00010"], digits=3))
  } else {
    message("  hsa00010: NOT enriched in GSEA")
  }
  
  # Dotplot
  gsea_kegg_sig <- gsea_kegg_df |> filter(p.adjust < 0.05)
  if (nrow(gsea_kegg_sig) > 0) {
    dp_df <- head(gsea_kegg_sig[order(gsea_kegg_sig$p.adjust), ], 20)
    dp_df$Description <- factor(dp_df$Description, levels = rev(dp_df$Description))
    p_gk <- ggplot(dp_df, aes(NES, Description, size = setSize, color = p.adjust)) +
      geom_point() + scale_color_gradient(low = "#8A2BE2", high = "#FFE135", trans = "log10") +
      labs(title = "GSEA KEGG — Significant Gene Sets", x = "NES") +
      theme_minimal(12) + theme(plot.background = element_rect(fill = "white", color = NA))
    ggsave(file.path(repo_root, "results", "figures", "GSEA_KEGG_dotplot.png"), p_gk, width = 10, height = 7, dpi = 300)
  }
}

# ── GSEA Reactome ──
gsea_reactome <- gsePathway(
  geneList = ranked_list,
  organism = "human",
  pvalueCutoff = 1
)

if (!is.null(gsea_reactome) && nrow(gsea_reactome) > 0) {
  gsea_react_df <- as_tibble(gsea_reactome@result)
  rio::export(gsea_react_df, file.path(repo_root, "results", "enrichment", "GSEA_Reactome.csv"))
  message("GSEA Reactome: ", nrow(gsea_react_df), " gene sets")
  
  # Check glycolysis-related
  glyc_re <- gsea_react_df[grepl("Glycolysis|Glucose|Carbohydrate|Pyruvate", gsea_react_df$Description, ignore.case = TRUE), ]
  if (nrow(glyc_re) > 0) {
    message("  Glycolysis-related Reactome terms found: ", nrow(glyc_re))
  } else {
    message("  No glycolysis-related terms significant in GSEA Reactome")
  }
  
  gsea_react_sig <- gsea_react_df |> filter(p.adjust < 0.05)
  if (nrow(gsea_react_sig) > 0) {
    dp_df <- head(gsea_react_sig[order(gsea_react_sig$p.adjust), ], 20)
    dp_df$Description <- factor(dp_df$Description, levels = rev(dp_df$Description))
    p_gr <- ggplot(dp_df, aes(NES, Description, size = setSize, color = p.adjust)) +
      geom_point() + scale_color_gradient(low = "#8A2BE2", high = "#FFE135", trans = "log10") +
      labs(title = "GSEA Reactome — Significant Gene Sets", x = "NES") +
      theme_minimal(12) + theme(plot.background = element_rect(fill = "white", color = NA))
    ggsave(file.path(repo_root, "results", "figures", "GSEA_Reactome_dotplot.png"), p_gr, width = 10, height = 7, dpi = 300)
  }
}

# ── Compare ORA vs GSEA ──
ora_kegg_file <- file.path(repo_root, "results", "enrichment", "KEGG_ORA_Up.csv")
if (file.exists(ora_kegg_file)) {
  ora_up <- rio::import(ora_kegg_file)
  comparison <- tibble(
    method = c("ORA_Up", "ORA_Down", "GSEA"),
    hsa00010_present = c(
      "hsa00010" %in% ora_up$ID,
      if(file.exists(file.path(repo_root,"results","enrichment","KEGG_ORA_Down.csv"))) "hsa00010" %in% rio::import(file.path(repo_root,"results","enrichment","KEGG_ORA_Down.csv"))$ID else FALSE,
      if(!is.null(gsea_kegg) && nrow(gsea_kegg)>0) "hsa00010" %in% gsea_kegg@result$ID else FALSE
    )
  )
  rio::export(comparison, file.path(repo_root, "results", "tables", "ora_vs_gsea_comparison.csv"))
}

message("\n✓ GSEA complete.")
