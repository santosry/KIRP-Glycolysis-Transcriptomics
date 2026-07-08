# 11_integrative_analysis.R
# INTEGRAÇÃO: Expressão Diferencial × ORA × PPI

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(rio)
  library(ggplot2)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
deg_file   <- file.path(repo_root, "results", "differential_expression", "DEG_global.csv")
cent_file  <- file.path(repo_root, "results", "tables", "ppi_centrality.csv")

dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)

# ── Load data ──
deg_all <- rio::import(deg_file)
cent <- if (file.exists(cent_file)) rio::import(cent_file) else NULL

# ── Build integrative table ──
# For DEGs that are also in PPI
if (!is.null(cent)) {
  integrative <- deg_all |>
    filter(gene_id %in% cent$gene_symbol) |>
    left_join(cent |> select(gene_symbol, degree, betweenness, closeness, eigenvector, community),
              by = c("gene_id" = "gene_symbol")) |>
    select(gene_id, logFC, adj.P.Val, regulation, AveExpr,
           degree, betweenness, closeness, eigenvector, community) |>
    arrange(desc(abs(logFC)))
} else {
  integrative <- deg_all |> filter(regulation %in% c("Up", "Down")) |>
    select(gene_id, logFC, adj.P.Val, regulation, AveExpr) |>
    arrange(desc(abs(logFC)))
}

# ── Add hsa00010 membership ──
kegg_file <- file.path(repo_root, "data", "metadata", "Hsa_genes.csv")
if (file.exists(kegg_file)) {
  hsa_genes <- rio::import(kegg_file)$gene_symbol
  integrative$in_hsa00010 <- integrative$gene_id %in% hsa_genes
}

# ── Add KEGG enrichment ──
kegg_up_file <- file.path(repo_root, "results", "enrichment", "KEGG_ORA_Up.csv")
kegg_down_file <- file.path(repo_root, "results", "enrichment", "KEGG_ORA_Down.csv")

if (file.exists(kegg_up_file)) {
  kegg_up <- rio::import(kegg_up_file)
  # For each gene, find which enriched KEGG pathways it belongs to
  if ("geneID" %in% colnames(kegg_up)) {
    integrative$kegg_up_pathways <- sapply(integrative$gene_id, function(g) {
      pw <- kegg_up$Description[grepl(g, kegg_up$geneID)]
      if (length(pw) == 0) return("") else return(paste(unique(pw), collapse = "; "))
    })
  }
}

# ── Priority score: simple convergence ──
# Number of evidence layers supporting a gene
integrative$evidence_layers <- 0
integrative$evidence_layers <- integrative$evidence_layers + 1  # DE
if (!is.null(cent)) {
  integrative$evidence_layers <- integrative$evidence_layers +
    (!is.na(integrative$degree))  # in PPI
  integrative$evidence_layers <- integrative$evidence_layers +
    (integrative$degree >= median(cent$degree[cent$gene_symbol %in% integrative$gene_id], na.rm = TRUE))  # high degree
}
integrative$evidence_layers <- integrative$evidence_layers +
  as.integer(integrative$in_hsa00010 & !is.na(integrative$in_hsa00010))  # in hsa00010

integrative <- integrative |> arrange(desc(evidence_layers), desc(abs(logFC)))

rio::export(integrative, file.path(repo_root, "results", "tables", "integrative_gene_table.csv"))

# ── Summary ──
message("\nTop genes by evidence convergence:")
top_int <- head(integrative, 10)
for (i in 1:nrow(top_int)) {
  r <- top_int[i, ]
  extras <- c()
  if (isTRUE(r$in_hsa00010)) extras <- c(extras, "hsa00010")
  if (!is.null(cent) && !is.na(r$degree)) extras <- c(extras, paste0("deg=", r$degree))
  message(sprintf("  %2d. %-10s logFC=% 6.2f  FDR=%.1e  %s  layers=%d  %s",
                  i, r$gene_id, r$logFC, r$adj.P.Val, r$regulation,
                  r$evidence_layers, paste(extras, collapse = ", ")))
}

message("\n✓ Integrative analysis complete. ", nrow(integrative), " genes in table.")
