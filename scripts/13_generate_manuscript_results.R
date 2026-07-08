# 13_generate_manuscript_results.R
# Gera manuscript_results.csv a partir dos outputs do pipeline

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
dir.create(file.path(repo_root, "results", "manuscript"), recursive = TRUE, showWarnings = FALSE)

read_if <- function(path) { if (file.exists(path)) rio::import(path) else NULL }

# ── Sample counts ──
meta <- readRDS(file.path(repo_root, "data", "processed", "metadata.rds"))
n_kirp <- sum(meta$condition == "KIRP")
n_gtx  <- sum(meta$condition == "Normal_GTEx")
n_tcga_norm <- 0
tcga_file <- file.path(repo_root, "data", "processed", "metadata_tcga_only.rds")
if (file.exists(tcga_file)) {
  tcga_meta <- readRDS(tcga_file)
  n_tcga_norm <- sum(tcga_meta$proposed_condition == "TCGA_Normal")
}

# ── DEG counts ──
deg_df <- read_if(file.path(repo_root, "results", "differential_expression", "DEG_global.csv"))
n_up <- sum(deg_df$regulation == "Up"); n_down <- sum(deg_df$regulation == "Down")

# ── PCA ──
pca_df <- read_if(file.path(repo_root, "results", "tables", "pca_variance.csv"))
pc1 <- if(!is.null(pca_df)) pca_df$variance_pct[1] else NA
pc2 <- if(!is.null(pca_df)) pca_df$variance_pct[2] else NA

# ── hsa00010 ──
hsa_audit <- read_if(file.path(repo_root, "results", "tables", "hsa00010_audit.csv"))
n_hsa_tested <- if(!is.null(hsa_audit)) sum(hsa_audit$tested, na.rm = TRUE) else NA
n_hsa_up     <- if(!is.null(hsa_audit)) sum(hsa_audit$regulation == "Up", na.rm = TRUE) else NA
n_hsa_down   <- if(!is.null(hsa_audit)) sum(hsa_audit$regulation == "Down", na.rm = TRUE) else NA

# ── ORA ──
kegg_up <- read_if(file.path(repo_root, "results", "enrichment", "KEGG_ORA_Up.csv"))
kegg_down <- read_if(file.path(repo_root, "results", "enrichment", "KEGG_ORA_Down.csv"))
top_kegg_up <- if(!is.null(kegg_up)) paste(head(kegg_up$Description[kegg_up$p.adjust < 0.05], 3), collapse = "; ") else ""
top_kegg_down <- if(!is.null(kegg_down)) paste(head(kegg_down$Description[kegg_down$p.adjust < 0.05], 3), collapse = "; ") else ""

react_up <- read_if(file.path(repo_root, "results", "enrichment", "Reactome_ORA_Up.csv"))
react_down <- read_if(file.path(repo_root, "results", "enrichment", "Reactome_ORA_Down.csv"))
top_react_up <- if(!is.null(react_up)) paste(head(react_up$Description[react_up$p.adjust < 0.05], 3), collapse = "; ") else ""
top_react_down <- if(!is.null(react_down)) paste(head(react_down$Description[react_down$p.adjust < 0.05], 3), collapse = "; ") else ""

# hsa00010 in ORA
hsa_kegg_up <- !is.null(kegg_up) && "hsa00010" %in% kegg_up$ID
hsa_kegg_down <- !is.null(kegg_down) && "hsa00010" %in% kegg_down$ID

# ── GSEA ──
gsea_kegg <- read_if(file.path(repo_root, "results", "enrichment", "GSEA_KEGG.csv"))
hsa_gsea <- !is.null(gsea_kegg) && "hsa00010" %in% gsea_kegg$ID

# ── STRING ──
net_sum <- read_if(file.path(repo_root, "results", "tables", "ppi_network_summary.csv"))
n_string_input  <- if(!is.null(net_sum)) net_sum$value[net_sum$metric == "n_DEGs_input"] else NA
n_string_mapped <- if(!is.null(net_sum)) net_sum$value[net_sum$metric == "n_mapped"] else NA
n_string_giant  <- if(!is.null(net_sum)) net_sum$value[net_sum$metric == "n_giant_component"] else NA
n_comm_louvain  <- if(!is.null(net_sum)) net_sum$value[net_sum$metric == "n_communities_Louvain"] else NA
mod_louvain     <- if(!is.null(net_sum)) net_sum$value[net_sum$metric == "modularity_Louvain"] else NA
ari_val         <- if(!is.null(net_sum)) net_sum$value[net_sum$metric == "ARI_Louvain_Leiden"] else NA

# ── Compile ──
mr <- tibble(
  key = c(
    "n_samples_KIRP", "n_samples_GTEx_Normal", "n_samples_TCGA_Normal",
    "n_genes_tested", "n_DEG_up", "n_DEG_down", "n_DEG_total",
    "PC1_variance_pct", "PC2_variance_pct",
    "n_hsa00010_in_matrix", "n_hsa00010_tested", "n_hsa00010_up", "n_hsa00010_down",
    "top_KEGG_Up_pathways", "top_KEGG_Down_pathways",
    "top_Reactome_Up_pathways", "top_Reactome_Down_pathways",
    "hsa00010_in_KEGG_ORA_Up", "hsa00010_in_KEGG_ORA_Down", "hsa00010_in_GSEA",
    "n_STRING_input", "n_STRING_mapped", "n_STRING_giant_component",
    "n_Louvain_communities", "modularity_Louvain", "ARI_Louvain_Leiden"
  ),
  value = c(
    n_kirp, n_gtx, n_tcga_norm,
    nrow(deg_df), n_up, n_down, n_up + n_down,
    pc1, pc2,
    if(!is.null(hsa_audit)) sum(hsa_audit$in_matrix, na.rm=TRUE) else NA,
    n_hsa_tested, n_hsa_up, n_hsa_down,
    top_kegg_up, top_kegg_down, top_react_up, top_react_down,
    hsa_kegg_up, hsa_kegg_down, hsa_gsea,
    n_string_input, n_string_mapped, n_string_giant,
    n_comm_louvain, mod_louvain, ari_val
  )
)
rio::export(mr, file.path(repo_root, "results", "manuscript", "manuscript_results.csv"))
message("✓ manuscript_results.csv generated with ", nrow(mr), " entries.")
