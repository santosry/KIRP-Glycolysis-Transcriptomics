# 10_string_network.R — STRING com pesos corrigidos (distance ≠ strength)

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
message("DEGs for STRING: ", length(all_deg), " (", length(up_genes), " Up, ", length(down_genes), " Down)")

if (length(all_deg) < 2) stop("Too few DEGs for PPI.")

node_annot <- tibble(gene_symbol = all_deg, regulation = ifelse(gene_symbol %in% up_genes, "Up", "Down"))

# STRING
options(timeout = 5000000)
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 0, input_directory = "")
mapped <- string_db$map(as.data.frame(node_annot), "gene_symbol", removeUnmappedRows = TRUE) |> filter(!is.na(STRING_id))
message("Mapped: ", nrow(mapped), " / ", length(all_deg))

unmapped <- setdiff(all_deg, mapped$gene_symbol)
if (length(unmapped) > 0) message("Unmapped: ", paste(unmapped, collapse = ", "))

score_cutoff <- 700
ppi <- string_db$get_interactions(mapped$STRING_id) |>
  filter(from %in% mapped$STRING_id, to %in% mapped$STRING_id, combined_score >= score_cutoff)
message("Edges (score ≥ ", score_cutoff, "): ", nrow(ppi))

rio::export(mapped, file.path(repo_root, "results", "ppi", "STRING_mapping.csv"))
rio::export(ppi, file.path(repo_root, "results", "ppi", "STRING_edges.csv"))

# Build graph
g <- graph_from_data_frame(
  d = ppi |> transmute(from, to, weight = combined_score),
  directed = FALSE,
  vertices = mapped |> transmute(name = STRING_id, gene_symbol, regulation)
)

# Largest component
comp <- components(g); giant <- which.max(comp$csize)
g_cc <- induced_subgraph(g, vids = V(g)[comp$membership == giant])
n_isolated <- vcount(g) - vcount(g_cc)
message("Largest component: ", vcount(g_cc), " nodes, ", ecount(g_cc), " edges | Isolated/removed: ", n_isolated)

# ═══════════════════════════════════════
# CORRECTED CENTRALITIES
# combined_score / 1000 = edge_strength (0-1)
# edge_distance = 1 / edge_strength (for weighted betweenness/closeness)
# ═══════════════════════════════════════
edge_strength <- E(g_cc)$weight / 1000
edge_distance <- 1 / edge_strength

V(g_cc)$degree              <- degree(g_cc)
V(g_cc)$strength            <- strength(g_cc, weights = edge_strength)
V(g_cc)$betweenness_unw     <- betweenness(g_cc, normalized = TRUE)
V(g_cc)$closeness_unw       <- closeness(g_cc, normalized = TRUE)
V(g_cc)$betweenness_w       <- betweenness(g_cc, weights = edge_distance, normalized = TRUE)
V(g_cc)$closeness_w         <- closeness(g_cc, weights = edge_distance, normalized = TRUE)
V(g_cc)$eigenvector_w       <- eigen_centrality(g_cc, weights = edge_strength)$vector

# ── Communities: Louvain + Leiden + ARI ──
set.seed(1)
comm_l <- cluster_louvain(g_cc)
comm_ld <- tryCatch(cluster_leiden(g_cc, objective_function = "modularity"), error = function(e) NULL)

mod_l <- modularity(comm_l)
mod_ld <- if(!is.null(comm_ld)) modularity(comm_ld) else NA

V(g_cc)$community_louvain <- as.character(membership(comm_l))
if (!is.null(comm_ld)) V(g_cc)$community_leiden <- as.character(membership(comm_ld))

# ARI between Louvain and Leiden
ari_val <- NA
if (!is.null(comm_ld)) {
  ari_val <- aricode::ARI(membership(comm_l), membership(comm_ld))
  message("ARI (Louvain vs Leiden): ", round(ari_val, 4))
}

# ── Centrality table ──
cent <- tibble(
  gene_symbol       = V(g_cc)$gene_symbol,
  regulation        = V(g_cc)$regulation,
  community_louvain = V(g_cc)$community_louvain,
  degree            = V(g_cc)$degree,
  strength          = round(V(g_cc)$strength, 3),
  betweenness_unw   = round(V(g_cc)$betweenness_unw, 4),
  closeness_unw     = round(V(g_cc)$closeness_unw, 4),
  betweenness_w     = round(V(g_cc)$betweenness_w, 4),
  closeness_w       = round(V(g_cc)$closeness_w, 4),
  eigenvector_w     = round(V(g_cc)$eigenvector_w, 4)
) |> arrange(desc(degree))
rio::export(cent, file.path(repo_root, "results", "tables", "ppi_centrality.csv"))

# ── Community summary ──
comm_sum <- cent |> group_by(community_louvain) |>
  summarise(n=n(), n_up=sum(regulation=="Up"), n_down=sum(regulation=="Down"),
            genes=paste(gene_symbol, collapse=", "), .groups="drop")
rio::export(comm_sum, file.path(repo_root, "results", "tables", "ppi_communities.csv"))

# ── Network visualization ──
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

# Community-colored
n_comm <- length(unique(V(g_cc)$community_louvain))
comm_cols <- colorRampPalette(brewer.pal(min(8, max(3, n_comm)), "Set1"))(n_comm)
set.seed(1)
p_comm <- ggraph(g_cc, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.15, color = "grey60") +
  scale_edge_width(range = c(0.2, 1.5), guide = "none") +
  geom_node_point(aes(color = community_louvain), size = 4, alpha = 0.9) +
  scale_color_manual(values = comm_cols) +
  geom_node_text(aes(label = gene_symbol), repel = TRUE, size = 3, max.overlaps = 50) +
  labs(color = "Community") + theme_void(14) +
  theme(plot.background = element_rect(fill = "#FFF8DC", color = NA))
ggsave(file.path(repo_root, "results", "figures", "PPI_communities.png"), p_comm, width = 10, height = 8, dpi = 300)

# ── Network summary ──
net_sum <- tibble(
  metric = c("n_DEGs_input", "n_mapped", "n_unmapped", "n_giant_component", "n_edges",
             "n_isolated_or_small", "n_communities_Louvain", "modularity_Louvain",
             "n_communities_Leiden", "modularity_Leiden", "ARI_Louvain_Leiden", "score_cutoff"),
  value = c(length(all_deg), nrow(mapped), length(unmapped), vcount(g_cc), ecount(g_cc),
            n_isolated, n_comm, round(mod_l, 4),
            if(!is.null(comm_ld)) length(unique(membership(comm_ld))) else NA,
            if(!is.null(comm_ld)) round(mod_ld, 4) else NA,
            round(ari_val, 4), score_cutoff)
)
rio::export(net_sum, file.path(repo_root, "results", "tables", "ppi_network_summary.csv"))

message("\n✓ STRING network complete. Communities: ", n_comm, " | Modularity: ", round(mod_l, 4))
