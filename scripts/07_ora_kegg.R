# 07_ora_kegg.R — ORA KEGG com mapeamento SYMBOL→ENTREZ unificado

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)  # provides AnnotationDbi::select for OrgDb queries
  library(enrichplot)
  library(dplyr)
  library(tibble)
  library(rio)
  library(ggplot2)
})
# org.Hs.eg.db masks dplyr::select — use dplyr::select for data frames

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
# UNIFIED SYMBOL → ENTREZ MAPPING
# ═══════════════════════════════════════
all_keys <- keys(org.Hs.eg.db, keytype = "SYMBOL")
valid_symbols <- intersect(universe_genes, all_keys)

map_raw <- AnnotationDbi::select(org.Hs.eg.db, keys = valid_symbols, columns = "ENTREZID", keytype = "SYMBOL")
map_raw <- map_raw[!is.na(map_raw$ENTREZID), ]

# Audit mappings
map_counts <- map_raw |> group_by(SYMBOL) |> summarise(n_entrez = n(), entrez_ids = paste(unique(ENTREZID), collapse = ";"), .groups = "drop")
map_summary <- map_counts |>
  mutate(mapping_status = case_when(n_entrez == 1 ~ "1:1", n_entrez > 1 ~ "1:n", TRUE ~ "unmapped"))

# Resolve 1:n: keep first Entrez ID (documented rule)
map_1to1 <- map_raw |> filter(SYMBOL %in% map_summary$SYMBOL[map_summary$mapping_status == "1:1"])
map_1ton <- map_raw |> filter(SYMBOL %in% map_summary$SYMBOL[map_summary$mapping_status == "1:n"]) |>
  group_by(SYMBOL) |> slice(1) |> ungroup()
gene_map <- bind_rows(map_1to1, map_1ton) |> dplyr::select(gene_symbol = SYMBOL, ENTREZID)
rio::export(gene_map, file.path(repo_root, "results", "tables", "gene_id_mapping.csv"))

unmapped_symbols <- setdiff(universe_genes, gene_map$gene_symbol)
message(sprintf("Gene mapping: %d 1:1 + %d 1:n = %d mapped | %d unmapped",
                nrow(map_1to1), nrow(map_1ton), nrow(gene_map), length(unmapped_symbols)))

# Universes for ORA
universe_entrez <- unique(gene_map$ENTREZID)
message("ORA universe: ", length(universe_entrez), " Entrez IDs")

# ═══════════════════════════════════════
# ORA KEGG: UP
# ═══════════════════════════════════════
up_symbols <- deg_all$gene_id[deg_all$regulation == "Up"]
up_entrez <- gene_map$ENTREZID[gene_map$gene_symbol %in% up_symbols]
message("Up DEGs mapped: ", length(up_entrez), " / ", length(up_symbols))

kegg_up <- NULL; kegg_up_sig <- NULL
if (length(up_entrez) >= 5) {
  kegg_up <- enrichKEGG(gene = up_entrez, universe = universe_entrez, organism = "hsa",
                         pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
  if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
    kegg_up_df <- as_tibble(kegg_up@result)
    rio::export(kegg_up_df, file.path(repo_root, "results", "enrichment", "KEGG_ORA_Up.csv"))
    kegg_up_sig <- kegg_up_df |> filter(p.adjust < 0.05)
    message("KEGG Up: ", nrow(kegg_up_sig), " significant / ", nrow(kegg_up_df), " total")
  }
}

# ORA KEGG: DOWN
down_symbols <- deg_all$gene_id[deg_all$regulation == "Down"]
down_entrez <- gene_map$ENTREZID[gene_map$gene_symbol %in% down_symbols]
message("Down DEGs mapped: ", length(down_entrez), " / ", length(down_symbols))

kegg_down <- NULL; kegg_down_sig <- NULL
if (length(down_entrez) >= 5) {
  kegg_down <- enrichKEGG(gene = down_entrez, universe = universe_entrez, organism = "hsa",
                           pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
  if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
    kegg_down_df <- as_tibble(kegg_down@result)
    rio::export(kegg_down_df, file.path(repo_root, "results", "enrichment", "KEGG_ORA_Down.csv"))
    kegg_down_sig <- kegg_down_df |> filter(p.adjust < 0.05)
    message("KEGG Down: ", nrow(kegg_down_sig), " significant / ", nrow(kegg_down_df), " total")
  }
}

# Dotplots
plot_ora <- function(df, dir_label) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df <- head(df[order(df$p.adjust), ], 20)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  ggplot(df, aes(Count, Description, size = Count, color = p.adjust)) +
    geom_point() + scale_color_gradient(low = "#8A2BE2", high = "#FFE135", trans = "log10", name = "FDR") +
    labs(title = paste0("KEGG ORA — ", dir_label), x = "Gene Count") +
    theme_minimal(12) + theme(plot.background = element_rect(fill = "#FFF8DC", color = NA))
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
  for (i in 1:nrow(ora_df)) {
    entrez_ids <- strsplit(ora_df$geneID[i], "/")[[1]]
    symbols <- gene_map$gene_symbol[match(entrez_ids, gene_map$ENTREZID)]
    symbols <- symbols[!is.na(symbols)]
    if (length(symbols) > 0) {
      rows[[length(rows) + 1]] <- tibble(
        database = db, direction = dir_label,
        pathway_id = ora_df$ID[i], pathway_name = ora_df$Description[i],
        pathway_fdr = ora_df$p.adjust[i], ENTREZID = entrez_ids,
        gene_symbol = symbols
      )
    }
  }
  if (length(rows) > 0) return(bind_rows(rows)) else return(NULL)
}

pw_kegg_up <- build_pw_table(kegg_up_df, "KEGG", "Up")
pw_kegg_down <- build_pw_table(kegg_down_df, "KEGG", "Down")
pw_all <- bind_rows(pw_kegg_up, pw_kegg_down)
if (!is.null(pw_all) && nrow(pw_all) > 0) {
  rio::export(pw_all, file.path(repo_root, "results", "tables", "pathway_gene_membership_KEGG.csv"))
}

message("\n✓ KEGG ORA complete. Universe: ", length(universe_entrez), " | Mapped: ", nrow(gene_map))
