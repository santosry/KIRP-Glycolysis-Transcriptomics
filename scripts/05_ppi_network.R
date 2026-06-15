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
deg_file <- file.path(repo_root, "results", "differential_expression", "DEG_KIRP_vs_Normal.csv")

if (!file.exists(deg_file)) {
  stop("Run scripts/03_differential_expression.R before generating the PPI network.")
}

dir.create(file.path(repo_root, "results", "ppi"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)

deg <- rio::import(deg_file)
up_genes <- deg |> filter(regulation == "Up") |> pull(gene_symbol) |> unique()
down_genes <- deg |> filter(regulation == "Down") |> pull(gene_symbol) |> unique()
deg_genes <- sort(unique(c(up_genes, down_genes)))

if (length(deg_genes) < 2) {
  stop("Too few DEGs for a meaningful PPI network.")
}

node_annot <- tibble(
  gene_symbol = deg_genes,
  regulation = case_when(
    gene_symbol %in% up_genes ~ "Up",
    gene_symbol %in% down_genes ~ "Down",
    TRUE ~ "NS"
  )
)

options(timeout = 5000000)
string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,
  score_threshold = 0,
  input_directory = ""
)

mapped <- string_db$map(as.data.frame(node_annot), "gene_symbol", removeUnmappedRows = TRUE)
mapped <- mapped |> filter(!is.na(STRING_id))
if (nrow(mapped) < 2) {
  stop("STRING mapping returned fewer than two proteins.")
}

mapped_ids <- mapped$STRING_id
ppi_raw <- string_db$get_interactions(mapped_ids)
score_cutoff <- 700

ppi <- ppi_raw |>
  filter(from %in% mapped_ids, to %in% mapped_ids) |>
  filter(combined_score >= score_cutoff)

if (nrow(ppi) == 0) {
  stop("No PPI edges after high-confidence filtering.")
}

rio::export(mapped, file.path(repo_root, "results", "ppi", "STRING_mapping.csv"))
rio::export(ppi, file.path(repo_root, "results", "ppi", "STRING_edges_high_confidence.csv"))

g <- graph_from_data_frame(
  d = ppi |> transmute(from, to, weight = combined_score),
  directed = FALSE,
  vertices = mapped |> transmute(name = STRING_id, gene_symbol = gene_symbol, regulation = regulation)
)

components_info <- components(g)
giant <- which.max(components_info$csize)
g_cc <- induced_subgraph(g, vids = V(g)[components_info$membership == giant])

set.seed(1)
p_ppi <- ggraph(g_cc, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.25) +
  scale_edge_width(range = c(0.2, 2.2), guide = "none") +
  geom_node_point(aes(color = regulation), size = 4) +
  scale_color_manual(values = c(Up = "#1F5BFF", Down = "#8A2BE2", NS = "grey70")) +
  geom_node_text(aes(label = gene_symbol), repel = TRUE, size = 3.5) +
  labs(color = "Regulation") +
  theme_void(base_size = 14)

ggsave(
  filename = file.path(repo_root, "results", "figures", "PPI_network.png"),
  plot = p_ppi,
  width = 9,
  height = 7,
  dpi = 300
)

network_summary <- tibble(
  metric = c("nodes_full", "edges_full", "nodes_largest_component", "edges_largest_component", "score_cutoff"),
  value = c(vcount(g), ecount(g), vcount(g_cc), ecount(g_cc), score_cutoff)
)
rio::export(network_summary, file.path(repo_root, "results", "tables", "ppi_network_summary.csv"))

message("Saved results/figures/PPI_network.png")
