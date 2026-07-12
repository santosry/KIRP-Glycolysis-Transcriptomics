# 07_ora_kegg.R — ORA KEGG com mapeamento SYMBOL→ENTREZ unificado
# (base R version to avoid dplyr / AnnotationDbi conflicts)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(rio)
  library(ggplot2)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
deg_file <- file.path(repo_root, "results", "differential_expression", "DEG_global.csv")
universe_file <- file.path(repo_root, "results", "differential_expression", "gene_universe.csv")
if (!file.exists(deg_file)) stop("Run 05 first.")

dir.create(file.path(repo_root, "results", "enrichment"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

deg_all <- rio::import(deg_file)
universe_df <- rio::import(universe_file)
universe_genes <- unique(universe_df$gene_id)

# ═══════════════════════════════════════
# UNIFIED SYMBOL → ENTREZ MAPPING (base R)
# ═══════════════════════════════════════
all_keys <- keys(org.Hs.eg.db, keytype = "SYMBOL")
valid_symbols <- intersect(universe_genes, all_keys)

map_raw <- AnnotationDbi::select(org.Hs.eg.db, keys = valid_symbols, columns = "ENTREZID", keytype = "SYMBOL")
# Force to plain data.frame with plain character columns
map_raw <- data.frame(
  SYMBOL = as.character(map_raw$SYMBOL),
  ENTREZID = as.character(map_raw$ENTREZID),
  stringsAsFactors = FALSE
)
map_raw <- map_raw[!is.na(map_raw$ENTREZID) & map_raw$ENTREZID != "", ]

# Audit mappings (base R)
sym_tab <- table(map_raw$SYMBOL)
n_1to1 <- sum(sym_tab == 1)
n_1ton <- sum(sym_tab > 1)
message(sprintf("Gene mapping: %d 1:1 + %d 1:n = %d mapped", n_1to1, n_1ton, length(sym_tab)))

# Resolve 1:n: keep first Entrez ID per symbol
keep_idx <- !duplicated(map_raw$SYMBOL)
gene_map <- map_raw[keep_idx, ]
names(gene_map) <- c("gene_symbol", "ENTREZID")
rio::export(gene_map, file.path(repo_root, "results", "tables", "gene_id_mapping.csv"))

unmapped_symbols <- setdiff(universe_genes, gene_map$gene_symbol)
message(sprintf("%d unmapped symbols", length(unmapped_symbols)))

# Universes for ORA
universe_entrez <- unique(gene_map$ENTREZID)
message("ORA universe: ", length(universe_entrez), " Entrez IDs")

# ═══════════════════════════════════════
# ORA KEGG: UP
# ═══════════════════════════════════════
up_symbols <- deg_all$gene_id[deg_all$regulation == "Up"]
up_entrez <- gene_map$ENTREZID[gene_map$gene_symbol %in% up_symbols]
message("Up DEGs mapped: ", length(up_entrez), " / ", length(up_symbols))

kegg_up <- NULL; kegg_up_df <- NULL; kegg_up_sig <- NULL
if (length(up_entrez) >= 5) {
  kegg_up <- enrichKEGG(gene = up_entrez, universe = universe_entrez, organism = "hsa",
                         pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
  if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
    kegg_up_df <- as.data.frame(kegg_up@result, stringsAsFactors = FALSE)
    rio::export(kegg_up_df, file.path(repo_root, "results", "enrichment", "KEGG_ORA_Up.csv"))
    kegg_up_sig <- kegg_up_df[kegg_up_df$p.adjust < 0.05, , drop = FALSE]
    message("KEGG Up: ", nrow(kegg_up_sig), " significant / ", nrow(kegg_up_df), " total")
  }
}

# ORA KEGG: DOWN
down_symbols <- deg_all$gene_id[deg_all$regulation == "Down"]
down_entrez <- gene_map$ENTREZID[gene_map$gene_symbol %in% down_symbols]
message("Down DEGs mapped: ", length(down_entrez), " / ", length(down_symbols))

kegg_down <- NULL; kegg_down_df <- NULL; kegg_down_sig <- NULL
if (length(down_entrez) >= 5) {
  kegg_down <- enrichKEGG(gene = down_entrez, universe = universe_entrez, organism = "hsa",
                           pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
  if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
    kegg_down_df <- as.data.frame(kegg_down@result, stringsAsFactors = FALSE)
    rio::export(kegg_down_df, file.path(repo_root, "results", "enrichment", "KEGG_ORA_Down.csv"))
    kegg_down_sig <- kegg_down_df[kegg_down_df$p.adjust < 0.05, , drop = FALSE]
    message("KEGG Down: ", nrow(kegg_down_sig), " significant / ", nrow(kegg_down_df), " total")
  }
}

# Dotplots
plot_ora <- function(df, dir_label) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- head(df[order(df$p.adjust), ], 20)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  ggplot(df, aes_string("Count", "Description", size = "Count", color = "p.adjust")) +
    geom_point() + scale_color_gradient(low = "#8A2BE2", high = "#FFE135", trans = "log10", name = "FDR") +
    labs(title = paste0("KEGG ORA — ", dir_label), x = "Gene Count") +
    theme_minimal(12) + theme(plot.background = element_rect(fill = "white", color = NA))
}
if (!is.null(kegg_up_sig)) ggsave(file.path(repo_root, "results", "figures", "KEGG_dotplot_Up.png"), plot_ora(kegg_up_sig, "Up"), width = 10, height = 7, dpi = 300)
if (!is.null(kegg_down_sig)) ggsave(file.path(repo_root, "results", "figures", "KEGG_dotplot_Down.png"), plot_ora(kegg_down_sig, "Down"), width = 10, height = 7, dpi = 300)

# Check hsa00010
for (dir_label in c("Up", "Down")) {
  df <- if(dir_label == "Up") kegg_up_df else if(dir_label == "Down") kegg_down_df else NULL
  present <- !is.null(df) && "hsa00010" %in% df$ID
  fdr_val <- if(present) format(df$p.adjust[df$ID == "hsa00010"], digits = 3) else "N/A"
  message(sprintf("  hsa00010 in KEGG %s: %s (FDR=%s)", dir_label, present, fdr_val))
}

# Build pathway-gene membership table for integration
build_pw_table <- function(ora_df, db, dir_label) {
  if (is.null(ora_df) || nrow(ora_df) == 0 || !("geneID" %in% colnames(ora_df))) return(NULL)
  rows <- list()
  for (i in seq_len(nrow(ora_df))) {
    entrez_ids <- strsplit(ora_df$geneID[i], "/")[[1]]
    symbols <- gene_map$gene_symbol[match(entrez_ids, gene_map$ENTREZID)]
    symbols <- symbols[!is.na(symbols)]
    if (length(symbols) > 0) {
      rows[[length(rows) + 1]] <- data.frame(
        database = db, direction = dir_label,
        pathway_id = ora_df$ID[i], pathway_name = ora_df$Description[i],
        pathway_fdr = ora_df$p.adjust[i], ENTREZID = entrez_ids,
        gene_symbol = symbols,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) > 0) return(do.call(rbind, rows)) else return(NULL)
}

pw_kegg_up <- build_pw_table(kegg_up_df, "KEGG", "Up")
pw_kegg_down <- build_pw_table(kegg_down_df, "KEGG", "Down")
pw_list <- list(pw_kegg_up, pw_kegg_down)
pw_list <- pw_list[!sapply(pw_list, is.null)]
if (length(pw_list) > 0) {
  pw_all <- do.call(rbind, pw_list)
  rio::export(pw_all, file.path(repo_root, "results", "tables", "pathway_gene_membership_KEGG.csv"))
}

message("\n✓ KEGG ORA complete. Universe: ", length(universe_entrez), " | Mapped: ", nrow(gene_map))
