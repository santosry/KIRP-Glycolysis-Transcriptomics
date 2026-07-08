# 08_ora_reactome.R
# ORA REACTOME — separadamente para Up e Down

suppressPackageStartupMessages({
  library(ReactomePA)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tibble)
  library(rio)
  library(ggplot2)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
deg_file <- file.path(repo_root, "results", "differential_expression", "DEG_global.csv")

if (!file.exists(deg_file)) stop("Run 05_differential_expression_global.R first.")

dir.create(file.path(repo_root, "results", "enrichment"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)

deg_all <- rio::import(deg_file)

# ── Build Entrez mapping ──
all_keys <- keys(org.Hs.eg.db, keytype = "SYMBOL")
all_genes <- unique(deg_all$gene_id)
valid_genes <- intersect(all_genes, all_keys)

symbol2entrez <- AnnotationDbi::select(org.Hs.eg.db, keys = valid_genes,
                         columns = "ENTREZID", keytype = "SYMBOL")
symbol2entrez <- symbol2entrez[!is.na(symbol2entrez$ENTREZID), ]
symbol2entrez <- symbol2entrez[!duplicated(symbol2entrez$SYMBOL), ]

universe_entrez <- unique(symbol2entrez$ENTREZID)
message("Universe Entrez: ", length(universe_entrez))

# ── ORA Reactome: UP ──
up_entrez <- symbol2entrez$ENTREZID[symbol2entrez$SYMBOL %in% deg_all$gene_id[deg_all$regulation == "Up"]]
message("Up Entrez: ", length(up_entrez))

reactome_up <- NULL
if (length(up_entrez) >= 5) {
  reactome_up <- enrichPathway(
    gene = up_entrez,
    universe = universe_entrez,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  if (!is.null(reactome_up) && nrow(reactome_up) > 0) {
    reactome_up <- as_tibble(reactome_up@result)
    rio::export(reactome_up, file.path(repo_root, "results", "enrichment", "Reactome_ORA_Up.csv"))
    message("Reactome Up: ", nrow(reactome_up), " terms")
  }
}

# ── ORA Reactome: DOWN ──
down_entrez <- symbol2entrez$ENTREZID[symbol2entrez$SYMBOL %in% deg_all$gene_id[deg_all$regulation == "Down"]]
message("Down Entrez: ", length(down_entrez))

reactome_down <- NULL
if (length(down_entrez) >= 5) {
  reactome_down <- enrichPathway(
    gene = down_entrez,
    universe = universe_entrez,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  if (!is.null(reactome_down) && nrow(reactome_down) > 0) {
    reactome_down <- as_tibble(reactome_down@result)
    rio::export(reactome_down, file.path(repo_root, "results", "enrichment", "Reactome_ORA_Down.csv"))
    message("Reactome Down: ", nrow(reactome_down), " terms")
  }
}

# ── Dotplots ──
plot_reactome <- function(df, direction, max_terms = 20) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  df_sig <- df |> filter(p.adjust < 0.05)
  if (nrow(df_sig) == 0) { message("No significant Reactome terms for ", direction); return(NULL) }
  df_sig <- head(df_sig[order(df_sig$p.adjust), ], max_terms)
  df_sig$Description <- factor(df_sig$Description, levels = rev(df_sig$Description))
  
  ggplot(df_sig, aes(x = Count, y = Description, size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#8A2BE2", high = "#FFE135", trans = "log10", name = "FDR") +
    labs(title = paste0("Reactome ORA — ", direction), x = "Gene Count", y = "") +
    theme_minimal(base_size = 12) +
    theme(plot.background = element_rect(fill = "#FFF8DC", color = NA))
}

if (!is.null(reactome_up)) {
  p <- plot_reactome(reactome_up, "Upregulated")
  if (!is.null(p)) ggsave(file.path(repo_root, "results", "figures", "Reactome_dotplot_Up.png"),
                           p, width = 10, height = 7, dpi = 300)
}
if (!is.null(reactome_down)) {
  p <- plot_reactome(reactome_down, "Downregulated")
  if (!is.null(p)) ggsave(file.path(repo_root, "results", "figures", "Reactome_dotplot_Down.png"),
                           p, width = 10, height = 7, dpi = 300)
}

# ── Check glycolysis-related Reactome terms ──
glyc_terms <- c("Glycolysis", "Gluconeogenesis", "Glucose metabolism",
                "Metabolism of carbohydrates", "Pyruvate metabolism")
for (gt in glyc_terms) {
  in_up <- !is.null(reactome_up) && any(grepl(gt, reactome_up$Description, ignore.case = TRUE))
  in_down <- !is.null(reactome_down) && any(grepl(gt, reactome_down$Description, ignore.case = TRUE))
  message(sprintf("  %-35s Up:%s Down:%s", gt, in_up, in_down))
}

# ── Concordance KEGG vs Reactome ──
concordance <- list()
if (!is.null(kegg_file <- file.path(repo_root, "results", "enrichment", "KEGG_ORA_Up.csv"))) {
  if (file.exists(kegg_file)) concordance$kegg_up <- rio::import(kegg_file)
}
if (!is.null(kegg_file2 <- file.path(repo_root, "results", "enrichment", "KEGG_ORA_Down.csv"))) {
  if (file.exists(kegg_file2)) concordance$kegg_down <- rio::import(kegg_file2)
}
concordance$reactome_up <- reactome_up
concordance$reactome_down <- reactome_down

# Build cross-ref table focusing on metabolism
meta_terms <- c(
  "carbohydrate|glycolysis|glucose|pyruvate|TCA|citrate|carbon|metabolic|amino acid|lipid|fatty acid"
)

concordance_table <- tibble(
  database = character(), direction = character(), term_id = character(),
  description = character(), fdr = numeric(), n_genes = integer()
)

for (db_name in names(concordance)) {
  df <- concordance[[db_name]]
  if (is.null(df) || nrow(df) == 0) next
  parts <- strsplit(db_name, "_")[[1]]
  direction <- parts[length(parts)]
  db <- paste(parts[-length(parts)], collapse = "_")
  
  df_filt <- df |> filter(grepl(meta_terms, Description, ignore.case = TRUE) & p.adjust < 0.05)
  if (nrow(df_filt) > 0) {
    concordance_table <- rbind(concordance_table, tibble(
      database = db, direction = direction,
      term_id = df_filt$ID, description = df_filt$Description,
      fdr = df_filt$p.adjust, n_genes = as.integer(df_filt$Count)
    ))
  }
}

if (nrow(concordance_table) > 0) {
  rio::export(concordance_table, file.path(repo_root, "results", "tables", "kegg_reactome_concordance.csv"))
}

message("\n✓ Reactome ORA complete.")
