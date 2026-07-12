# 16b_ppi_3d_correlation.R
# PPI/Co-expression Network 3D baseado em correlação (Pearson) entre os 106 genes
# Usa matrizes de expressão por via (RDS) ou transcriptoma completo
# 100% offline, reprodutível
#
# Output: results/v3/figures/PPI_network_3D.html

suppressPackageStartupMessages({
  library(plotly)
  library(dplyr)
  library(rio)
  library(igraph)
  library(htmlwidgets)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
dir.create(file.path(repo_root, "results", "v3", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "v3", "tables"), recursive = TRUE, showWarnings = FALSE)

message("=== 3D Co-expression Network ===")

# ── Load pathway expression matrices (RDS: samples x genes) ──
message("  Loading pathway expression matrices...")

pw_expr <- list()
pw_labels <- c("hsa00010", "hsa00020", "hsa00030")
all_gene_lists <- list()

for (pw_id in pw_labels) {
  rds_file <- file.path(repo_root, "data", "processed", paste0("expression_", pw_id, ".rds"))
  if (file.exists(rds_file)) {
    mat <- readRDS(rds_file)
    
    # Determine orientation: genes should be columns (genes < samples)
    # If rows > cols, then samples x genes → transpose
    if (nrow(mat) > ncol(mat)) {
      # mat is samples x genes → transpose to genes x samples
      E_pw <- t(mat)
    } else {
      E_pw <- mat
    }
    
    pw_expr[[pw_id]] <- E_pw
    all_gene_lists[[pw_id]] <- rownames(E_pw)
    message(sprintf("    %s: %d genes x %d samples", pw_id, nrow(E_pw), ncol(E_pw)))
  }
}

if (length(pw_expr) == 0) {
  stop("No pathway expression matrices found. Run scripts/00c_three_pathway_analysis.R first.")
}

# Combine all unique genes
all_pw_genes <- unique(unlist(all_gene_lists))
n_total_genes <- length(all_pw_genes)
message(sprintf("  Total unique pathway genes: %d", n_total_genes))

# Find common samples across all pathway matrices
common_samples <- Reduce(intersect, lapply(pw_expr, colnames))
message(sprintf("  Common samples across pathways: %d", length(common_samples)))

# Build combined matrix (genes x samples)
E_combined <- matrix(NA, nrow = n_total_genes, ncol = length(common_samples))
rownames(E_combined) <- all_pw_genes
colnames(E_combined) <- common_samples

for (pw_id in pw_labels) {
  if (pw_id %in% names(pw_expr)) {
    mat <- pw_expr[[pw_id]]
    genes_in_pw <- intersect(rownames(mat), all_pw_genes)
    samps_in_pw <- intersect(colnames(mat), common_samples)
    E_combined[genes_in_pw, samps_in_pw] <- mat[genes_in_pw, samps_in_pw]
  }
}

storage.mode(E_combined) <- "numeric"
message(sprintf("  Combined matrix: %d genes x %d samples", nrow(E_combined), ncol(E_combined)))

# ── Load DEG results ──
deg_all <- do.call(rbind, lapply(pw_labels, function(pw_id) {
  f <- file.path(repo_root, "results", "v3", "tables", paste0("DEG_", pw_id, ".csv"))
  if (file.exists(f)) {
    d <- rio::import(f)
    d$pathway <- pw_id
    return(d)
  }
  return(NULL)
}))

deg_all <- deg_all[!duplicated(deg_all$gene_id), ]
deg_all$regulation <- "NS"
deg_all$regulation[deg_all$FDR_Paired < 0.05 & deg_all$logFC_Paired > 1] <- "Up"
deg_all$regulation[deg_all$FDR_Paired < 0.05 & deg_all$logFC_Paired < -1] <- "Down"

n_up   <- sum(deg_all$regulation == "Up", na.rm = TRUE)
n_down <- sum(deg_all$regulation == "Down", na.rm = TRUE)
message(sprintf("  DEGs: %d Up, %d Down, %d NS", n_up, n_down,
                sum(deg_all$regulation == "NS")))

# ── Correlation network ──
message("  Computing Pearson correlations...")
cor_mat <- cor(t(E_combined), method = "pearson", use = "pairwise.complete.obs")

cor_threshold <- 0.6
message(sprintf("  Threshold: |r| > %.1f", cor_threshold))

# Build edges
edges_list <- list()
for (i in 1:(nrow(cor_mat) - 1)) {
  for (j in (i + 1):nrow(cor_mat)) {
    r_ij <- cor_mat[i, j]
    if (!is.na(r_ij) && abs(r_ij) > cor_threshold) {
      edges_list[[length(edges_list) + 1]] <- data.frame(
        from = rownames(cor_mat)[i],
        to   = rownames(cor_mat)[j],
        weight = abs(r_ij),
        correlation = r_ij,
        stringsAsFactors = FALSE
      )
    }
  }
}
edges_df <- do.call(rbind, edges_list)

# Fallback threshold
if (is.null(edges_df) || nrow(edges_df) == 0) {
  cor_threshold <- 0.4
  message(sprintf("  No edges at 0.6, lowering to %.1f", cor_threshold))
  edges_list <- list()
  for (i in 1:(nrow(cor_mat) - 1)) {
    for (j in (i + 1):nrow(cor_mat)) {
      r_ij <- cor_mat[i, j]
      if (!is.na(r_ij) && abs(r_ij) > cor_threshold) {
        edges_list[[length(edges_list) + 1]] <- data.frame(
          from = rownames(cor_mat)[i],
          to   = rownames(cor_mat)[j],
          weight = abs(r_ij),
          correlation = r_ij,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  edges_df <- do.call(rbind, edges_list)
}

message(sprintf("  Edges: %d", nrow(edges_df)))

# ── Node attributes ──
node_df <- data.frame(
  name = all_pw_genes,
  gene_symbol = all_pw_genes,
  stringsAsFactors = FALSE
)
node_df$regulation <- deg_all$regulation[match(node_df$gene_symbol, deg_all$gene_id)]
node_df$regulation[is.na(node_df$regulation)] <- "NS"
node_df$logFC <- deg_all$logFC_Paired[match(node_df$gene_symbol, deg_all$gene_id)]
node_df$FDR   <- deg_all$FDR_Paired[match(node_df$gene_symbol, deg_all$gene_id)]

# ── Build graph ──
g <- graph_from_data_frame(
  d = edges_df[, c("from", "to")],
  directed = FALSE,
  vertices = node_df
)
E(g)$weight <- edges_df$weight

comp <- components(g)
giant_idx <- which.max(comp$csize)
g_cc <- induced_subgraph(g, vids = V(g)[comp$membership == giant_idx])
message(sprintf("  Largest component: %d nodes, %d edges", vcount(g_cc), ecount(g_cc)))

# ── 3D Layout ──
set.seed(42)
layout_3d <- layout_with_fr(g_cc, dim = 3, niter = 5000)

# Metrics
V(g_cc)$degree <- degree(g_cc)
V(g_cc)$betweenness <- betweenness(g_cc, normalized = TRUE)

node_out <- data.frame(
  gene_symbol = V(g_cc)$gene_symbol,
  regulation  = V(g_cc)$regulation,
  logFC       = round(V(g_cc)$logFC, 3),
  degree      = V(g_cc)$degree,
  betweenness = round(V(g_cc)$betweenness, 4),
  x = layout_3d[, 1],
  y = layout_3d[, 2],
  z = layout_3d[, 3],
  stringsAsFactors = FALSE
)

# Node sizes by degree
if (max(node_out$degree) > min(node_out$degree)) {
  node_out$size <- 5 + 12 * (node_out$degree - min(node_out$degree)) /
    (max(node_out$degree) - min(node_out$degree))
} else {
  node_out$size <- 8
}

# ── Edge traces ──
el <- as_edgelist(g_cc)
name_to_idx <- setNames(seq_len(nrow(node_out)), V(g_cc)$name)

edge_x <- edge_y <- edge_z <- c()
for (i in seq_len(nrow(el))) {
  fi <- name_to_idx[el[i, 1]]
  ti <- name_to_idx[el[i, 2]]
  if (!is.na(fi) && !is.na(ti)) {
    edge_x <- c(edge_x, node_out$x[fi], node_out$x[ti], NA)
    edge_y <- c(edge_y, node_out$y[fi], node_out$y[ti], NA)
    edge_z <- c(edge_z, node_out$z[fi], node_out$z[ti], NA)
  }
}

# ── Build figure ──
message("  Building 3D plotly...")

fig_ppi <- plot_ly() %>%
  add_trace(
    x = edge_x, y = edge_y, z = edge_z,
    type = "scatter3d", mode = "lines",
    line = list(color = "rgba(150,150,150,0.25)", width = 0.6),
    hoverinfo = "none", showlegend = FALSE
  )

for (reg in c("Up", "Down", "NS")) {
  nd <- node_out[node_out$regulation == reg, ]
  if (nrow(nd) == 0) next
  
  col <- c(Up = "#1F5BFF", Down = "#8A2BE2", NS = "grey70")[reg]
  lbl <- c(Up = "\u25B2 Up", Down = "\u25BC Down", NS = "NS")[reg]
  op  <- c(Up = 0.95, Down = 0.95, NS = 0.6)[reg]
  
  fig_ppi <- fig_ppi %>%
    add_trace(
      data = nd,
      x = ~x, y = ~y, z = ~z,
      type = "scatter3d", mode = "markers+text",
      marker = list(
        size = ~size, color = col, opacity = op,
        line = list(color = "rgba(0,0,0,0.4)", width = 0.8)
      ),
      text = ~gene_symbol,
      textposition = "top center",
      textfont = list(size = 10, color = col),
      name = lbl,
      hoverinfo = "text",
      hovertext = ~paste0(
        "<b>", gene_symbol, "</b><br>",
        "Regulation: ", regulation, "<br>",
        "log2FC: ", logFC, "<br>",
        "Degree: ", degree, "<br>",
        "Betweenness: ", betweenness
      )
    )
}

fig_ppi <- fig_ppi %>%
  layout(
    title = list(
      text = paste0(
        "<b>Gene Co-expression Network \u2014 Central Carbon Metabolism (3D)</b><br>",
        "<sup>Pearson |r| > ", cor_threshold,
        " | ", vcount(g_cc), " nodes | ", ecount(g_cc), " edges | ",
        "Fruchterman-Reingold 3D</sup>"
      ),
      font = list(size = 16)
    ),
    scene = list(
      xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      zaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
      camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))
    ),
    showlegend = TRUE,
    legend = list(title = list(text = "<b>Regulation</b>"), x = 0.8, y = 0.95),
    margin = list(l = 0, r = 0, b = 0, t = 50)
  ) %>%
  config(
    displayModeBar = TRUE,
    modeBarButtonsToRemove = c("sendDataToCloud", "autoScale2d",
                                "hoverClosestCartesian", "hoverCompareCartesian",
                                "toggleSpikelines"),
    displaylogo = FALSE, scrollZoom = TRUE, doubleClick = "reset"
  )

out_ppi <- file.path(repo_root, "results", "v3", "figures", "PPI_network_3D.html")
saveWidget(fig_ppi, file = out_ppi, selfcontained = FALSE,
           title = "3D Co-expression Network \u2014 Central Carbon Metabolism in KIRP")
message(sprintf("  Saved: %s", out_ppi))

# ── Save tables ──
rio::export(node_out, file.path(repo_root, "results", "v3", "tables", "ppi_3d_centrality.csv"))

net_sum <- data.frame(
  metric = c("genes_total", "nodes_giant", "edges",
             "cor_threshold", "method", "n_up", "n_down", "n_ns"),
  value = c(n_total_genes, vcount(g_cc), ecount(g_cc),
            cor_threshold, "pearson_abs",
            sum(node_out$regulation == "Up"),
            sum(node_out$regulation == "Down"),
            sum(node_out$regulation == "NS"))
)
rio::export(net_sum, file.path(repo_root, "results", "v3", "tables", "ppi_3d_summary.csv"))

message("=== 3D PPI Complete ===")
