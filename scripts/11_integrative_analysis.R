# 11_integrative_analysis.R
# Integração DE × ORA × PPI — sem score arbitrário, convergência descritiva

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
deg_file   <- file.path(repo_root, "results", "differential_expression", "DEG_global.csv")
cent_file  <- file.path(repo_root, "results", "tables", "ppi_centrality.csv")
pw_file    <- file.path(repo_root, "results", "tables", "pathway_gene_membership_KEGG.csv")
hsa_file   <- file.path(repo_root, "data", "metadata", "Hsa_genes.csv")

dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

deg_all <- rio::import(deg_file) |> filter(regulation %in% c("Up", "Down"))
cent <- if (file.exists(cent_file)) rio::import(cent_file) else NULL
pw_kegg <- if (file.exists(pw_file)) rio::import(pw_file) else NULL
hsa_genes <- if (file.exists(hsa_file)) rio::import(hsa_file)$gene_symbol else NULL

# ── Build per-gene table ──
int_table <- deg_all |>
  select(gene_id, logFC, adj.P.Val, regulation, AveExpr, t) |>
  rename(gene_symbol = gene_id, moderated_t = t) |>
  mutate(in_hsa00010 = gene_symbol %in% hsa_genes)

# Merge KEGG pathways
if (!is.null(pw_kegg)) {
  pw_summary <- pw_kegg |>
    group_by(gene_symbol, direction) |>
    summarise(
      KEGG_pathways = paste(unique(pathway_name), collapse = "; "),
      n_KEGG_pathways = n_distinct(pathway_id),
      .groups = "drop"
    )
  
  pw_up <- pw_summary |> filter(direction == "Up") |> select(-direction) |>
    rename(KEGG_pathways_Up = KEGG_pathways, n_KEGG_Up = n_KEGG_pathways)
  pw_down <- pw_summary |> filter(direction == "Down") |> select(-direction) |>
    rename(KEGG_pathways_Down = KEGG_pathways, n_KEGG_Down = n_KEGG_pathways)
  
  int_table <- int_table |>
    left_join(pw_up, by = "gene_symbol") |>
    left_join(pw_down, by = "gene_symbol")
  
  int_table$n_KEGG_Up[is.na(int_table$n_KEGG_Up)] <- 0
  int_table$n_KEGG_Down[is.na(int_table$n_KEGG_Down)] <- 0
}

# Merge Reactome pathways (if available)
react_file <- file.path(repo_root, "results", "tables", "pathway_gene_membership_Reactome.csv")
if (file.exists(react_file)) {
  pw_react <- rio::import(react_file)
  pw_r_sum <- pw_react |> group_by(gene_symbol, direction) |>
    summarise(Reactome_pathways = paste(unique(pathway_name), collapse = "; "),
              n_Reactome = n_distinct(pathway_id), .groups = "drop")
  
  r_up <- pw_r_sum |> filter(direction == "Up") |> select(-direction) |>
    rename(Reactome_pathways_Up = Reactome_pathways, n_Reactome_Up = n_Reactome)
  r_down <- pw_r_sum |> filter(direction == "Down") |> select(-direction) |>
    rename(Reactome_pathways_Down = Reactome_pathways, n_Reactome_Down = n_Reactome)
  
  int_table <- int_table |> left_join(r_up, by = "gene_symbol") |> left_join(r_down, by = "gene_symbol")
}

# Merge PPI centrality
if (!is.null(cent)) {
  int_table <- int_table |> left_join(cent |> select(gene_symbol, community_louvain,
    degree, strength, betweenness_unw, closeness_unw, betweenness_w, closeness_w, eigenvector_w),
    by = "gene_symbol")
}

# Sort by |logFC| (largest effect first), no priority score
int_table <- int_table |> arrange(desc(abs(logFC)))

rio::export(int_table, file.path(repo_root, "results", "tables", "integrative_gene_table.csv"))

# Descriptive summary
message("\nIntegrative table: ", nrow(int_table), " genes")
if (!is.null(cent)) {
  in_ppi <- sum(!is.na(int_table$community_louvain))
  message("  In PPI largest component: ", in_ppi)
}
if (!is.null(hsa_genes)) {
  in_hsa <- sum(int_table$in_hsa00010, na.rm = TRUE)
  message("  In hsa00010: ", in_hsa)
}

# Top 10 by |logFC| with their evidence
message("\nTop 10 genes by |logFC|:")
for (i in 1:min(10, nrow(int_table))) {
  r <- int_table[i, ]
  details <- c()
  if (isTRUE(r$in_hsa00010)) details <- c(details, "hsa00010")
  if (!is.null(cent) && !is.na(r$degree)) details <- c(details, paste0("deg=", r$degree))
  if (!is.null(r$n_KEGG_Up) && r$n_KEGG_Up > 0) details <- c(details, paste0("KEGGup=", r$n_KEGG_Up))
  message(sprintf("  %2d. %-10s logFC=% 7.2f FDR=%.1e %s [%s]",
                  i, r$gene_symbol, r$logFC, r$adj.P.Val, r$regulation, paste(details, collapse=", ")))
}

message("\n✓ Integrative analysis complete.")
