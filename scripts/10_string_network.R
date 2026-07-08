# 10_string_network.R
# REDE STRING RECONSTRUÍDA com comunidades (Leiden/Louvain) e centralidade

suppressPackageStartupMessages({
  library(STRINGdb)
  library(dplyr)
  library(tibble)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
deg_file <- file.path(repo_root, "results", "differential_expression", "DEG_global.csv")

if (!file.exists(deg_file)) stop("Run 05_differential_expression_global.R first.")

dir.create(file.path(repo_root, "results", "ppi"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

deg_all <- rio::import(deg_file)

# ── Input: DEGs (FDR < 0.05, |logFC| > 1) ──
deg_genes <- deg_all |> filter(regulation %in% c("Up", "Down"))
up_genes <- unique(deg_genes$gene_id[deg_genes$regulation == "Up"])
down_genes <- unique(deg_genes$gene_id[deg_genes$regulation == "Down"])
all_deg_symbols <- sort(unique(c(up_genes, down_genes)))

message("DEGs submitted to STRING: ", length(all_deg_symbols))
message("  Up: ", length(up_genes), " | Down: ", length(down_genes))

if (length(all_deg_symbols) < 2) stop("Too few DEGs for PPI network.")

node_annot <- tibble(
  gene_symbol = all_deg_symbols,
  regulation = ifelse(gene_symbol %in% up_genes, "Up", "Down")
)

# ── STRING mapping ──
options(timeout = 5000000)
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 0, input_directory = "")

mapped <- string_db$map(as.data.frame(node_annot), "gene_symbol", removeUnmappedRows = TRUE)
mapped <- mapped |> filter(!is.na(STRING_id))
message("Mapped to STRING: ", nrow(mapped), " / ", length(all_deg_symbols))

unmapped <- setdiff(all_deg_symbols, mapped$gene_symbol)
if (length(unmapped) > 0) message("Unmapped: ", paste(unmapped, collapse = ", "))

if (nrow(mapped) < 2) stop("Too few mapped proteins.")

# ── Get interactions ──
ppi_raw <- string_db$get_interactions(mapped$STRING_id)
score_cutoff <- 700

ppi <- ppi_raw |>
  filter(from %in% mapped$STRING_id, to %in% mapped$STRING_id) |>
  filter(combined_score >= score_cutoff)

message("Edges (score ≥ ", score_cutoff, "): ", nrow(ppi))

rio::export(mapped, file.path(repo_root, "results", "ppi", "STRING_mapping.csv"))
rio::export(ppi, file.path(repo_root, "results", "ppi", "STRING_edges.csv"))

# ── Build graph ──
g <- graph_from_data_frame(
  d = ppi |> transmute(from, to, weight = combined_score),
  directed = FALSE,
  vertices = mapped |> transmute(name = STRING_id, gene_symbol, regulation)
)

# Largest component
comp <- components(g)
giant <- which.max(comp$csize)
g_cc <- induced_subgraph(g, vids = V(g)[comp$membership == giant])
message("Largest component: ", vcount(g_cc), " nodes, ", ecount(g_cc), " edges")

# ── Community detection (Louvain) ──
set.seed(1)
comm_louvain <- cluster_louvain(g_cc)
V(g_cc)$community_louvain <- as.character(membership(comm_louvain))

# ── Community detection (Leiden) ──
set.seed(1)
comm_leiden <- cluster_leiden(g_cc, objective_function = "modularity")
V(g_cc)$community_leiden <- as.character(membership(comm_leiden))

n_comm <- length(unique(V(g_cc)$community_louvain))
mod_louvain <- round(modularity(comm_louvain), 4)
mod_leiden <- round(modularity(comm_leiden), 4)
message("Communities (Louvain): ", n_comm, " | modularity: ", mod_louvain)
message("Communities (Leiden): ", length(unique(V(g_cc)$community_leiden)), " | modularity: ", mod_leiden)

# ── Centrality ──
V(g_cc)$degree      <- degree(g_cc)
V(g_cc)$betweenness <- betweenness(g_cc, normalized = TRUE)
V(g_cc)$closeness   <- closeness(g_cc, normalized = TRUE)
V(g_cc)$eigenvector <- eigen_centrality(g_cc)$vector

centrality_table <- tibble(
  gene_symbol = V(g_cc)$gene_symbol,
  regulation   = V(g_cc)$regulation,
  community    = V(g_cc)$community_louvain,
  degree       = V(g_cc)$degree,
  betweenness  = round(V(g_cc)$betweenness, 4),
  closeness    = round(V(g_cc)$closeness, 4),
  eigenvector  = round(V(g_cc)$eigenvector, 4)
) |> arrange(desc(degree))

rio::export(centrality_table, file.path(repo_root, "results", "tables", "ppi_centrality.csv"))

# Top by degree
message("\nTop genes by degree:")
top_deg <- head(centrality_table, 5)
for (i in 1:nrow(top_deg)) {
  message(sprintf("  %s: degree=%d, betweenness=%.3f, eigenvector=%.3f",
                  top_deg$gene_symbol[i], top_deg$degree[i],
                  top_deg$betweenness[i], top_deg$eigenvector[i]))
}

# ── Community summary ──
comm_summary <- centrality_table |>
  group_by(community) |>
  summarise(
    n_genes = n(), n_up = sum(regulation == "Up"), n_down = sum(regulation == "Down"),
    genes = paste(gene_symbol, collapse = ", "), .groups = "drop"
  )
rio::export(comm_summary, file.path(repo_root, "results", "tables", "ppi_communities.csv"))
message("\nCommunity summary:")
print(as.data.frame(comm_summary))

# ── Visualization ──
set.seed(1)
p_net <- ggraph(g_cc, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.2, color = "grey50") +
  scale_edge_width(range = c(0.2, 2.2), guide = "none") +
  geom_node_point(aes(color = regulation, size = degree), alpha = 0.9) +
  scale_color_manual(values = c(Up = "#1F5BFF", Down = "#8A2BE2")) +
  scale_size(range = c(2, 8)) +
  geom_node_text(aes(label = gene_symbol), repel = TRUE, size = 3, max.overlaps = 50) +
  labs(color = "Regulation", size = "Degree") +
  theme_void(base_size = 14) +
  theme(plot.background = element_rect(fill = "#FFF8DC", color = NA))

ggsave(file.path(repo_root, "results", "figures", "PPI_network.png"), p_net, width = 10, height = 8, dpi = 300)

# ── Community-colored network ──
n_comm_plot <- min(n_comm, 8)
comm_colors <- RColorBrewer::brewer.pal(max(3, n_comm_plot), "Set1")
if (n_comm > length(comm_colors)) comm_colors <- colorRampPalette(comm_colors)(n_comm)

set.seed(1)
p_comm <- ggraph(g_cc, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.15, color = "grey60") +
  scale_edge_width(range = c(0.2, 1.5), guide = "none") +
  geom_node_point(aes(color = community_louvain), size = 4, alpha = 0.9) +
  scale_color_manual(values = comm_colors) +
  geom_node_text(aes(label = gene_symbol), repel = TRUE, size = 3, max.overlaps = 50) +
  labs(color = "Comunidade") +
  theme_void(base_size = 14) +
  theme(plot.background = element_rect(fill = "#FFF8DC", color = NA))

ggsave(file.path(repo_root, "results", "figures", "PPI_communities.png"), p_comm, width = 10, height = 8, dpi = 300)

# ── Network summary ──
net_summary <- tibble(
  metric = c("nodes_input", "nodes_mapped", "nodes_unmapped",
             "nodes_giant_component", "edges_giant_component",
             "n_communities_louvain", "modularity_louvain",
             "n_communities_leiden", "modularity_leiden",
             "score_cutoff"),
  value = c(length(all_deg_symbols), nrow(mapped), length(unmapped),
            vcount(g_cc), ecount(g_cc),
            n_comm, mod_louvain,
            length(unique(V(g_cc)$community_leiden)), mod_leiden,
            score_cutoff)
)
rio::export(net_summary, file.path(repo_root, "results", "tables", "ppi_network_summary.csv"))

message("\n✓ STRING network analysis complete.")
