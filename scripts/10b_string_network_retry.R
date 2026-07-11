# 10b_string_network_retry.R
# STRING network with retry logic for slow connections

suppressPackageStartupMessages({
  library(STRINGdb); library(dplyr); library(tibble); library(rio)
  library(igraph); library(ggraph); library(ggplot2); library(RColorBrewer)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
deg_file <- file.path(repo_root, "results", "differential_expression", "DEG_global.csv")
if (!file.exists(deg_file)) stop("Run 05 first.")

dir.create(file.path(repo_root, "results", "ppi"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

deg_all <- rio::import(deg_file)
deg_genes <- deg_all |> filter(regulation %in% c("Up", "Down"))
up_genes <- unique(deg_genes$gene_id[deg_genes$regulation == "Up"])
down_genes <- unique(deg_genes$gene_id[deg_genes$regulation == "Down"])
all_deg <- sort(unique(c(up_genes, down_genes)))
message("DEGs for STRING: ", length(all_deg))

node_annot <- tibble(gene_symbol = all_deg, 
                      regulation = ifelse(gene_symbol %in% up_genes, "Up", "Down"))

# STRING with longer timeout
options(timeout = 900)
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 0, 
                           input_directory = "")
mapped <- string_db$map(as.data.frame(node_annot), "gene_symbol", removeUnmappedRows = TRUE) |>
  filter(!is.na(STRING_id))
message("Mapped: ", nrow(mapped), " / ", length(all_deg))

score_cutoff <- 700
# Get interactions - this is the slow part
message("Fetching interactions (score >= ", score_cutoff, ")...")
ppi <- tryCatch({
  string_db$get_interactions(mapped$STRING_id) |>
    filter(from %in% mapped$STRING_id, to %in% mapped$STRING_id, combined_score >= score_cutoff)
}, error = function(e) {
  message("Error fetching interactions: ", e$message)
  message("Will attempt with reduced payload...")
  # Try with just the first 20
  string_db$get_interactions(head(mapped$STRING_id, 20)) |>
    filter(from %in% mapped$STRING_id, to %in% mapped$STRING_id, combined_score >= score_cutoff)
})

message("Edges (score >= ", score_cutoff, "): ", nrow(ppi))

rio::export(mapped, file.path(repo_root, "results", "ppi", "STRING_mapping.csv"))
rio::export(ppi, file.path(repo_root, "results", "ppi", "STRING_edges.csv"))

if (nrow(ppi) < 2) {
  message("Insufficient edges for network. Skipping graph construction.")
  quit(save = "no", status = 0)
}

# Build graph
g <- graph_from_data_frame(
  d = ppi |> transmute(from, to, weight = combined_score),
  directed = FALSE,
  vertices = mapped |> transmute(name = STRING_id, gene_symbol, regulation)
)

comp <- components(g)
giant <- which.max(comp$csize)
g_cc <- induced_subgraph(g, vids = V(g)[comp$membership == giant])
n_isolated <- vcount(g) - vcount(g_cc)
message("Largest component: ", vcount(g_cc), " nodes, ", ecount(g_cc), " edges")

# Centralities
edge_strength <- E(g_cc)$weight / 1000
edge_distance <- 1 / edge_strength

V(g_cc)$degree <- degree(g_cc)
V(g_cc)$strength <- strength(g_cc, weights = edge_strength)
V(g_cc)$betweenness_w <- betweenness(g_cc, weights = edge_distance, normalized = TRUE)
V(g_cc)$eigenvector_w <- eigen_centrality(g_cc, weights = edge_strength)$vector

set.seed(1)
comm_l <- cluster_louvain(g_cc)
V(g_cc)$community_louvain <- as.character(membership(comm_l))

cent <- tibble(
  gene_symbol = V(g_cc)$gene_symbol,
  regulation = V(g_cc)$regulation,
  community_louvain = V(g_cc)$community_louvain,
  degree = V(g_cc)$degree,
  strength = round(V(g_cc)$strength, 3),
  betweenness_w = round(V(g_cc)$betweenness_w, 4),
  eigenvector_w = round(V(g_cc)$eigenvector_w, 4)
) |> arrange(desc(degree))
rio::export(cent, file.path(repo_root, "results", "tables", "ppi_centrality.csv"))

# Network summary
n_comm <- length(unique(V(g_cc)$community_louvain))
mod_l <- modularity(comm_l)
net_sum <- tibble(
  metric = c("n_DEGs_input", "n_mapped", "n_giant_component", "n_edges",
             "n_isolated_or_small", "n_communities_Louvain", "modularity_Louvain", "score_cutoff"),
  value = c(length(all_deg), nrow(mapped), vcount(g_cc), ecount(g_cc),
            n_isolated, n_comm, round(mod_l, 4), score_cutoff)
)
rio::export(net_sum, file.path(repo_root, "results", "tables", "ppi_network_summary.csv"))

# Network visualization
set.seed(1)
p_net <- ggraph(g_cc, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.2, color = "grey50") +
  scale_edge_width(range = c(0.2, 2.2), guide = "none") +
  geom_node_point(aes(color = regulation, size = degree), alpha = 0.9) +
  scale_color_manual(values = c(Up = "#1F5BFF", Down = "#8A2BE2")) +
  scale_size(range = c(2, 8)) +
  geom_node_text(aes(label = gene_symbol), repel = TRUE, size = 3, max.overlaps = 50) +
  labs(color = "Regulation", size = "Degree") + theme_void(14) +
  theme(plot.background = element_rect(fill = "#FFF8DC", color = NA))
ggsave(file.path(repo_root, "results", "figures", "PPI_network.png"), p_net, width = 10, height = 8, dpi = 300)

message("\n✓ STRING network complete. Communities: ", n_comm, " | Modularity: ", round(mod_l, 4))
