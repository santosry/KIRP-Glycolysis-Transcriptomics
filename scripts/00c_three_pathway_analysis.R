# 00c_three_pathway_analysis.R
# Complete analysis of all 3 central carbon metabolism pathways in KIRP
# For each pathway: DE, Volcano, Enrichment (ORA+GSEA), PPI network

suppressPackageStartupMessages({
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(rio)
  library(igraph)
  library(ggraph)
  library(RColorBrewer)
  library(STRINGdb)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)

# ── Load full data ──
expr_full <- readRDS(file.path(repo_root, "data", "processed", "expression_matrix_full.rds"))
meta <- readRDS(file.path(repo_root, "data", "processed", "metadata_full.rds"))

# Main analysis: KIRP vs GTEx Normal (same as original design)
main_idx <- meta$condition %in% c("KIRP", "Normal_GTEx")
meta_main <- meta[main_idx, ]
expr_main <- expr_full[main_idx, ]

message(sprintf("Main analysis: %d samples (%d KIRP + %d GTEx)", 
                nrow(meta_main), sum(meta_main$condition == "KIRP"), 
                sum(meta_main$condition == "Normal_GTEx")))

# ── Pathway definitions ──
pw <- rio::import(file.path(repo_root, "results", "tables", "pathway_gene_membership.csv"))
all_genes <- colnames(expr_full)

pathways <- list(
  hsa00010 = list(
    name = "Glycolysis / Gluconeogenesis",
    short = "Glycolysis",
    genes = intersect(pw$gene[pw$in_hsa00010], all_genes)
  ),
  hsa00030 = list(
    name = "Pentose Phosphate Pathway",
    short = "PPP",
    genes = intersect(pw$gene[pw$in_hsa00030], all_genes)
  ),
  hsa00020 = list(
    name = "Citrate Cycle (TCA)",
    short = "TCA",
    genes = intersect(pw$gene[pw$in_hsa00020], all_genes)
  )
)

# ── Output directories ──
for (d in c("results/pathways", "results/pathways/figures", "results/pathways/tables")) {
  dir.create(file.path(repo_root, d), recursive = TRUE, showWarnings = FALSE)
}

# ═════════════════════════════════════════════════════════════
# PER-PATHWAY ANALYSIS
# ═════════════════════════════════════════════════════════════

all_degs <- list()
all_ora_up <- list()
all_ora_down <- list()
all_gsea <- list()

for (pw_id in names(pathways)) {
  pw_info <- pathways[[pw_id]]
  pw_genes <- pw_info$genes
  n_genes <- length(pw_genes)
  
  message(sprintf("\n══════════ %s (%s): %d genes ══════════", pw_id, pw_info$short, n_genes))
  
  if (n_genes < 3) {
    message("  Too few genes, skipping.")
    next
  }
  
  # ── 1. Expression matrix for this pathway ──
  E <- t(expr_main[, pw_genes, drop = FALSE])
  storage.mode(E) <- "numeric"
  
  # ── 2. Differential Expression ──
  meta_pw <- meta_main
  meta_pw$condition <- factor(meta_pw$condition, levels = c("Normal_GTEx", "KIRP"))
  design <- model.matrix(~ 0 + condition, data = meta_pw)
  colnames(design) <- levels(meta_pw$condition)
  
  fit <- lmFit(E, design)
  fit <- eBayes(fit, robust = TRUE, trend = TRUE)
  
  contrast <- makeContrasts(KIRP - Normal_GTEx, levels = design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2, robust = TRUE, trend = TRUE)
  
  deg <- topTable(fit2, number = Inf, adjust.method = "BH")
  deg$gene_id <- rownames(deg)
  deg$regulation <- ifelse(deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1,
                           ifelse(deg$logFC > 0, "Up", "Down"), "NS")
  deg$pathway <- pw_id
  
  n_up <- sum(deg$regulation == "Up")
  n_down <- sum(deg$regulation == "Down")
  message(sprintf("  DEGs: %d Up | %d Down | %d NS", n_up, n_down, sum(deg$regulation == "NS")))
  
  rio::export(deg, file.path(repo_root, "results", "pathways", paste0("DEG_", pw_id, ".csv")))
  all_degs[[pw_id]] <- deg
  
  # ── 3. Volcano Plot ──
  deg_plot <- deg
  deg_plot$label <- ifelse(deg_plot$regulation != "NS", deg_plot$gene_id, "")
  
  p_volc <- ggplot(deg_plot, aes(logFC, -log10(adj.P.Val), color = regulation, label = label)) +
    geom_point(alpha = 0.7, size = 2.5) +
    scale_color_manual(values = c(Up = "#1F5BFF", Down = "#8A2BE2", NS = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 30, box.padding = 0.3) +
    labs(title = paste0(pw_info$short, " — KIRP vs Normal"),
         subtitle = sprintf("%d genes | %d Up | %d Down (FDR<0.05, |logFC|>1)", n_genes, n_up, n_down),
         x = "log2(Fold Change)", y = "-log10(FDR)") +
    theme_minimal(14) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(file.path(repo_root, "results", "pathways", "figures", paste0("Volcano_", pw_id, ".png")),
         p_volc, width = 10, height = 8, dpi = 300)
  
  # ── 4. Heatmap of DEGs ──
  deg_genes <- deg$gene_id[deg$regulation != "NS"]
  if (length(deg_genes) >= 3) {
    expr_sub <- expr_main[, deg_genes, drop = FALSE]
    gene_vars <- apply(expr_sub, 2, var, na.rm = TRUE)
    ok_genes <- deg_genes[gene_vars > 0.001 & is.finite(gene_vars)]
    if (length(ok_genes) >= 3) {
      deg_genes <- ok_genes
      heat_data <- t(scale(t(expr_main[, deg_genes, drop = FALSE])))
      heat_data[!is.finite(heat_data)] <- 0
    } else {
      deg_genes <- character(0)
    }
  }
  if (length(deg_genes) >= 3) {
    tryCatch({
    # Clip for visualization
    heat_data[heat_data > 3] <- 3
    heat_data[heat_data < -3] <- -3
    
    ann_col <- data.frame(
      Condition = meta_main$condition,
      row.names = rownames(meta_main)
    )
    ann_colors <- list(Condition = c(KIRP = "#8A2BE2", Normal_GTEx = "#1F5BFF"))
    
    png(file.path(repo_root, "results", "pathways", "figures", paste0("Heatmap_", pw_id, ".png")),
        width = 1200, height = 1000, res = 150)
    pheatmap::pheatmap(heat_data,
                        annotation_col = ann_col,
                        annotation_colors = ann_colors,
                        show_rownames = length(deg_genes) <= 50,
                        show_colnames = FALSE,
                        main = paste0(pw_info$short, " \U2014 ", length(deg_genes), " DEGs"),
                        color = colorRampPalette(c("#1F5BFF", "white", "#8A2BE2"))(100),
                        border_color = NA)
    dev.off()
    message("  Heatmap saved")
    }, error = function(e) message("  Heatmap error: ", e$message))
  }
  
  # ── 5. ORA Enrichment (KEGG) ──
  # Map to Entrez
  valid_genes <- intersect(pw_genes, keys(org.Hs.eg.db, keytype = "SYMBOL"))
  if (length(valid_genes) > 0) {
    map <- AnnotationDbi::select(org.Hs.eg.db, keys = valid_genes, columns = "ENTREZID", keytype = "SYMBOL")
    map <- data.frame(SYMBOL = as.character(map$SYMBOL), ENTREZID = as.character(map$ENTREZID), stringsAsFactors = FALSE)
    map <- map[!is.na(map$ENTREZID) & map$ENTREZID != "", ]
    map <- map[!duplicated(map$SYMBOL), ]
    
    universe_entrez <- unique(map$ENTREZID)
    
    # ORA Up
    up_symbols <- deg$gene_id[deg$regulation == "Up"]
    up_entrez <- map$ENTREZID[map$SYMBOL %in% up_symbols]
    
    if (length(up_entrez) >= 5) {
      ora_up <- enrichKEGG(gene = up_entrez, universe = universe_entrez, organism = "hsa",
                            pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
      if (!is.null(ora_up) && nrow(ora_up) > 0) {
        ora_up_df <- as.data.frame(ora_up@result)
        rio::export(ora_up_df, file.path(repo_root, "results", "pathways", paste0("ORA_Up_", pw_id, ".csv")))
        all_ora_up[[pw_id]] <- ora_up_df
        sig_up <- sum(ora_up_df$p.adjust < 0.05)
        message(sprintf("  ORA Up: %d sig / %d total", sig_up, nrow(ora_up_df)))
      }
    }
    
    # ORA Down
    down_symbols <- deg$gene_id[deg$regulation == "Down"]
    down_entrez <- map$ENTREZID[map$SYMBOL %in% down_symbols]
    
    if (length(down_entrez) >= 5) {
      ora_down <- enrichKEGG(gene = down_entrez, universe = universe_entrez, organism = "hsa",
                              pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
      if (!is.null(ora_down) && nrow(ora_down) > 0) {
        ora_down_df <- as.data.frame(ora_down@result)
        rio::export(ora_down_df, file.path(repo_root, "results", "pathways", paste0("ORA_Down_", pw_id, ".csv")))
        all_ora_down[[pw_id]] <- ora_down_df
        sig_down <- sum(ora_down_df$p.adjust < 0.05)
        message(sprintf("  ORA Down: %d sig / %d total", sig_down, nrow(ora_down_df)))
        
        # Dotplot for significant down
        if (sig_down > 0) {
          df <- ora_down_df[ora_down_df$p.adjust < 0.05, ]
          df <- head(df[order(df$p.adjust), ], 20)
          if (nrow(df) > 0) {
            df$Description <- factor(df$Description, levels = rev(df$Description))
            p_dot <- ggplot(df, aes(Count, Description, size = Count, color = p.adjust)) +
              geom_point() + 
              scale_color_gradient(low = "#8A2BE2", high = "#FFE135", trans = "log10", name = "FDR") +
              labs(title = paste0("KEGG ORA Down — ", pw_info$short), x = "Gene Count") +
              theme_minimal(12) +
              theme(plot.background = element_rect(fill = "white", color = NA))
            ggsave(file.path(repo_root, "results", "pathways", "figures", paste0("Dotplot_", pw_id, ".png")),
                   p_dot, width = 10, height = 7, dpi = 300)
          }
        }
      }
    }
  }
  
  # ── 6. GSEA ──
  t_stats <- deg$t
  names(t_stats) <- deg$gene_id
  t_stats <- sort(t_stats, decreasing = TRUE)
  
  if (length(t_stats) >= 5) {
    gsea_kegg <- gseKEGG(geneList = t_stats, organism = "hsa", pvalueCutoff = 1, eps = 0)
    if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
      gsea_df <- as.data.frame(gsea_kegg@result)
      rio::export(gsea_df, file.path(repo_root, "results", "pathways", paste0("GSEA_", pw_id, ".csv")))
      all_gsea[[pw_id]] <- gsea_df
      
      # Check if the pathway itself is enriched
      self_enriched <- pw_id %in% gsea_df$ID
      self_nes <- if(self_enriched) round(gsea_df$NES[gsea_df$ID == pw_id], 2) else NA
      self_fdr <- if(self_enriched) format(gsea_df$p.adjust[gsea_df$ID == pw_id], digits=2) else "NS"
      message(sprintf("  GSEA: %s in results: %s (NES=%s, FDR=%s)", pw_id, self_enriched, self_nes, self_fdr))
    }
  }
  
  # ── 7. PPI Network (STRING) - optional, skip on failure ──
  deg_only <- deg[deg$regulation != "NS", ]
  if (nrow(deg_only) >= 5) {
    tryCatch({
      message(sprintf("  Building PPI for %d DEGs...", nrow(deg_only)))
      
      node_annot <- data.frame(
        gene_symbol = deg_only$gene_id,
        regulation = deg_only$regulation,
        stringsAsFactors = FALSE
      )
      
      options(timeout = 120)
      string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 0, input_directory = "")
      mapped <- string_db$map(node_annot, "gene_symbol", removeUnmappedRows = TRUE)
      mapped <- mapped[!is.na(mapped$STRING_id), ]
      message(sprintf("    Mapped: %d / %d", nrow(mapped), nrow(deg_only)))
      
      if (nrow(mapped) >= 3) {
        score_cutoff <- 400
        ppi <- string_db$get_interactions(mapped$STRING_id)
        ppi <- ppi[ppi$from %in% mapped$STRING_id & ppi$to %in% mapped$STRING_id & ppi$combined_score >= score_cutoff, ]
        message(sprintf("    Edges (score>=%d): %d", score_cutoff, nrow(ppi)))
        
        if (nrow(ppi) >= 2) {
          g <- graph_from_data_frame(
            d = data.frame(from = ppi$from, to = ppi$to, weight = ppi$combined_score),
            directed = FALSE,
            vertices = data.frame(name = mapped$STRING_id, gene_symbol = mapped$gene_symbol, 
                                  regulation = mapped$regulation, stringsAsFactors = FALSE)
          )
          
          comp <- components(g)
          giant <- which.max(comp$csize)
          g_cc <- induced_subgraph(g, vids = V(g)[comp$membership == giant])
          V(g_cc)$degree <- degree(g_cc)
          
          set.seed(1)
          p_ppi <- ggraph(g_cc, layout = "fr") +
            geom_edge_link(aes(width = weight), alpha = 0.3, color = "grey50") +
            scale_edge_width(range = c(0.2, 2), guide = "none") +
            geom_node_point(aes(color = regulation, size = degree), alpha = 0.9) +
            scale_color_manual(values = c(Up = "#1F5BFF", Down = "#8A2BE2")) +
            scale_size(range = c(2, 7)) +
            geom_node_text(aes(label = gene_symbol), repel = TRUE, size = 3, max.overlaps = 50) +
            labs(color = "Regulation", size = "Degree",
                 title = paste0("PPI Network \U2014 ", pw_info$short)) +
            theme_void(14) +
            theme(plot.background = element_rect(fill = "white", color = NA))
          
          ggsave(file.path(repo_root, "results", "pathways", "figures", paste0("PPI_", pw_id, ".png")),
                 p_ppi, width = 10, height = 8, dpi = 300)
          
          ppi_out <- data.frame(
            gene_symbol = V(g_cc)$gene_symbol,
            regulation = V(g_cc)$regulation,
            degree = V(g_cc)$degree,
            stringsAsFactors = FALSE
          )
          rio::export(ppi_out, file.path(repo_root, "results", "pathways", paste0("PPI_", pw_id, ".csv")))
          message(sprintf("    PPI: %d nodes, %d edges", vcount(g_cc), ecount(g_cc)))
        }
      }
    }, error = function(e) {
      message(sprintf("    PPI skipped: %s", e$message))
    })
  }
}

# ═════════════════════════════════════════════════════════════
# Cross-pathway comparison
# ═════════════════════════════════════════════════════════════

message("\n══════════ CROSS-PATHWAY COMPARISON ══════════")

# Combine all DEGs
all_deg_combined <- do.call(rbind, all_degs)
all_deg_combined$gene_id <- rownames(all_deg_combined)
rownames(all_deg_combined) <- NULL

# Summary table
pw_summary_table <- data.frame(
  pathway = names(pathways),
  name = sapply(pathways, `[[`, "name"),
  genes_tested = sapply(pathways, function(p) length(p$genes)),
  deg_up = sapply(all_degs, function(d) sum(d$regulation == "Up")),
  deg_down = sapply(all_degs, function(d) sum(d$regulation == "Down")),
  deg_total = sapply(all_degs, function(d) sum(d$regulation != "NS")),
  stringsAsFactors = FALSE
)

# Add top DEG
for (i in seq_len(nrow(pw_summary_table))) {
  pw_id <- pw_summary_table$pathway[i]
  if (!is.null(all_degs[[pw_id]]) && nrow(all_degs[[pw_id]]) > 0) {
    d <- all_degs[[pw_id]]
    d_sig <- d[d$regulation != "NS", ]
    if (nrow(d_sig) > 0) {
      top_idx <- which.max(abs(d_sig$logFC))
      pw_summary_table$top_gene[i] <- d_sig$gene_id[top_idx]
      pw_summary_table$top_logFC[i] <- round(d_sig$logFC[top_idx], 2)
      pw_summary_table$top_regulation[i] <- d_sig$regulation[top_idx]
    }
  }
}

rio::export(pw_summary_table, file.path(repo_root, "results", "pathways", "pathway_summary.csv"))
print(pw_summary_table)

message("\nThree-pathway analysis COMPLETE!")
message(sprintf("Output directory: %s", file.path(repo_root, "results", "pathways")))
