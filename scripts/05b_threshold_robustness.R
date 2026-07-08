# 05b_threshold_robustness.R
# Executa ORA KEGG e Reactome nos 3 thresholds de DEG

suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(ReactomePA)
  library(dplyr); library(tibble); library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
deg_file <- file.path(repo_root, "results", "differential_expression", "DEG_global.csv")
map_file <- file.path(repo_root, "results", "tables", "gene_id_mapping.csv")
if (!file.exists(deg_file) || !file.exists(map_file)) stop("Run 05 and 07 first.")

dir.create(file.path(repo_root, "results", "enrichment"), recursive = TRUE, showWarnings = FALSE)

deg_all <- rio::import(deg_file)
gene_map <- rio::import(map_file)
universe_entrez <- unique(gene_map$ENTREZID)

thresholds <- list(
  "LFC0"  = 0,
  "LFC05" = 0.5,
  "LFC1"  = 1.0
)

robustness_rows <- list()

for (th_name in names(thresholds)) {
  lfc <- thresholds[[th_name]]
  
  # Classify at this threshold
  deg_all$regulation_th <- "NS"
  deg_all$regulation_th[deg_all$adj.P.Val < 0.05 & deg_all$logFC > lfc] <- "Up"
  deg_all$regulation_th[deg_all$adj.P.Val < 0.05 & deg_all$logFC < -lfc] <- "Down"
  
  for (direction in c("Up", "Down")) {
    genes <- deg_all$gene_id[deg_all$regulation_th == direction]
    entrez <- gene_map$ENTREZID[gene_map$gene_symbol %in% genes]
    
    if (length(entrez) < 5) next
    
    # KEGG
    kegg_res <- enrichKEGG(entrez, universe = universe_entrez, organism = "hsa",
                            pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
    if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
      kegg_df <- as_tibble(kegg_res@result)
      for (i in 1:nrow(kegg_df)) {
        robustness_rows[[length(robustness_rows) + 1]] <- tibble(
          database = "KEGG", direction = direction, threshold = th_name,
          pathway_id = kegg_df$ID[i], pathway_name = kegg_df$Description[i],
          FDR = kegg_df$p.adjust[i], significant = kegg_df$p.adjust[i] < 0.05,
          gene_count = as.integer(kegg_df$Count[i]),
          genes = kegg_df$geneID[i]
        )
      }
    }
    
    # Reactome
    react_res <- enrichPathway(entrez, universe = universe_entrez, organism = "human",
                                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
    if (!is.null(react_res) && nrow(react_res) > 0) {
      react_df <- as_tibble(react_res@result)
      for (i in 1:nrow(react_df)) {
        robustness_rows[[length(robustness_rows) + 1]] <- tibble(
          database = "Reactome", direction = direction, threshold = th_name,
          pathway_id = react_df$ID[i], pathway_name = react_df$Description[i],
          FDR = react_df$p.adjust[i], significant = react_df$p.adjust[i] < 0.05,
          gene_count = as.integer(react_df$Count[i]),
          genes = react_df$geneID[i]
        )
      }
    }
  }
}

if (length(robustness_rows) > 0) {
  robustness_df <- bind_rows(robustness_rows)
  rio::export(robustness_df, file.path(repo_root, "results", "tables", "pathway_threshold_robustness.csv"))
  
  # Summary: how many pathways are stable across thresholds?
  stable_summary <- robustness_df |>
    group_by(database, direction, pathway_id, pathway_name) |>
    summarise(
      significant_LFC0 = any(threshold == "LFC0" & significant),
      significant_LFC05 = any(threshold == "LFC05" & significant),
      significant_LFC1 = any(threshold == "LFC1" & significant),
      n_thresholds_sig = sum(significant),
      .groups = "drop"
    )
  rio::export(stable_summary, file.path(repo_root, "results", "tables", "pathway_stability_summary.csv"))
  message("Threshold robustness: ", nrow(robustness_df), " pathway-threshold combinations evaluated.")
}

# ── JACCARD KEGG vs REACTOME ──
# Compare gene sets of significant pathways
jaccard_rows <- list()
kegg_files <- c(
  file.path(repo_root, "results", "enrichment", "KEGG_ORA_Up.csv"),
  file.path(repo_root, "results", "enrichment", "KEGG_ORA_Down.csv")
)
react_files <- c(
  file.path(repo_root, "results", "enrichment", "Reactome_ORA_Up.csv"),
  file.path(repo_root, "results", "enrichment", "Reactome_ORA_Down.csv")
)

for (dir_idx in 1:2) {
  kf <- kegg_files[dir_idx]; rf <- react_files[dir_idx]
  if (!file.exists(kf) || !file.exists(rf)) next
  
  kegg_df <- rio::import(kf) |> filter(p.adjust < 0.05)
  react_df <- rio::import(rf) |> filter(p.adjust < 0.05)
  
  if (nrow(kegg_df) == 0 || nrow(react_df) == 0) next
  
  # Parse gene sets using the pathway_gene_membership table
  pw_file <- file.path(repo_root, "results", "tables", "pathway_gene_membership_KEGG.csv")
  if (!file.exists(pw_file)) next
  pw_all <- rio::import(pw_file)
  
  direction_label <- c("Up", "Down")[dir_idx]
  pw_dir <- pw_all |> filter(direction == direction_label)
  
  # Get unique pathways
  kegg_pathways <- unique(pw_dir$pathway_id)
  
  for (kp in kegg_pathways) {
    kegg_genes <- unique(pw_dir$gene_symbol[pw_dir$pathway_id == kp])
    
    for (rp in unique(pw_dir$pathway_id)) {
      if (kp == rp) next
      react_genes <- unique(pw_dir$gene_symbol[pw_dir$pathway_id == rp])
      
      intersection <- length(intersect(kegg_genes, react_genes))
      union_size <- length(union(kegg_genes, react_genes))
      jaccard <- if (union_size > 0) intersection / union_size else 0
      
      if (jaccard > 0.1) {  # Only report meaningful overlaps
        jaccard_rows[[length(jaccard_rows) + 1]] <- tibble(
          direction = direction_label,
          pathway_A = kp,
          pathway_B = rp,
          jaccard = round(jaccard, 4),
          n_intersection = intersection,
          n_union = union_size,
          genes_intersection = paste(intersect(kegg_genes, react_genes), collapse = ", ")
        )
      }
    }
  }
}

if (length(jaccard_rows) > 0) {
  jaccard_df <- bind_rows(jaccard_rows) |> arrange(desc(jaccard))
  rio::export(jaccard_df, file.path(repo_root, "results", "tables", "kegg_reactome_jaccard.csv"))
  message("Jaccard KEGG-Reactome: ", nrow(jaccard_df), " pathway pairs with Jaccard > 0.1")
}

message("\n✓ Threshold robustness + Jaccard analysis complete.")
