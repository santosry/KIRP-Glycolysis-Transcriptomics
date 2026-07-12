# 16_3d_visualizations.R
# Gera visualizações 3D interativas em HTML:
#   1. PPI network 3D (STRING interactome, todos os 106 genes + DEGs)
#   2. Volcano plots 3D (logFC × -log10(FDR) × AveExpr) para cada via
#
# Dependências: plotly, dplyr, rio, igraph, htmlwidgets, STRINGdb, scales
# Output: results/v3/figures/PPI_network_3D.html
#         results/v3/figures/Volcano_hsa00010_3D.html
#         results/v3/figures/Volcano_hsa00020_3D.html
#         results/v3/figures/Volcano_hsa00030_3D.html

suppressPackageStartupMessages({
  library(plotly)
  library(dplyr)
  library(rio)
  library(igraph)
  library(htmlwidgets)
  library(STRINGdb)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
dir.create(file.path(repo_root, "results", "v3", "figures"), recursive = TRUE, showWarnings = FALSE)

# ═══════════════════════════════════════════════════════════════════
# PART 1: 3D Volcano Plots (para cada via metabólica)
# ═══════════════════════════════════════════════════════════════════
message("=== 3D Volcano Plots ===")

pw_labels <- c(
  hsa00010 = "Glycolysis / Gluconeogenesis",
  hsa00020 = "Citrate Cycle (TCA)",
  hsa00030 = "Pentose Phosphate Pathway"
)

for (pw_id in names(pw_labels)) {
  deg_file <- file.path(repo_root, "results", "v3", "tables", paste0("DEG_", pw_id, ".csv"))
  if (!file.exists(deg_file)) {
    message(sprintf("  SKIP %s: file not found (%s)", pw_id, deg_file))
    next
  }
  
  deg <- rio::import(deg_file)
  
  # Classificação
  deg$regulation <- "NS"
  deg$regulation[deg$FDR_Paired < 0.05 & deg$logFC_Paired > 1] <- "Up"
  deg$regulation[deg$FDR_Paired < 0.05 & deg$logFC_Paired < -1] <- "Down"
  
  # -log10(FDR)
  deg$neg_log10_FDR <- -log10(pmax(deg$FDR_Paired, 1e-300, na.rm = TRUE))
  
  # Label
  deg$label <- ifelse(deg$regulation != "NS", deg$gene_id, "")
  
  n_up   <- sum(deg$regulation == "Up", na.rm = TRUE)
  n_down <- sum(deg$regulation == "Down", na.rm = TRUE)
  n_total <- nrow(deg)
  
  # ── Build Plotly ──
  fig <- plot_ly(
    data = deg,
    x = ~logFC_Paired,
    y = ~neg_log10_FDR,
    z = ~AveExpr,
    color = ~regulation,
    colors = c(Up = "#1F5BFF", Down = "#8A2BE2", NS = "grey70"),
    type = "scatter3d",
    mode = "markers",
    marker = list(
      size = 6,
      opacity = 0.85,
      line = list(color = "rgba(0,0,0,0.3)", width = 0.5),
      sizemode = "diameter"
    ),
    hoverinfo = "text",
    hovertext = ~paste0(
      "<b>", gene_id, "</b><br>",
      "log2FC: ", round(logFC_Paired, 3), "<br>",
      "FDR: ", formatC(FDR_Paired, format = "e", digits = 2), "<br>",
      "CI 95%: [", round(CI.L_Paired, 2), ", ", round(CI.R_Paired, 2), "]<br>",
      "AveExpr: ", round(AveExpr, 2), "<br>",
      "Regulation: ", regulation
    )
  ) %>%
    layout(
      title = list(
        text = paste0(
          "<b>", pw_labels[pw_id], "</b><br>",
          "<sup>3D Volcano — Paired Analysis (32 pairs) | ",
          n_total, " genes | ", n_up, "\u25B2 Up | ", n_down, " \u25BC Down</sup>"
        ),
        font = list(size = 16)
      ),
      scene = list(
        xaxis = list(title = "log<sub>2</sub>(Fold Change)", zeroline = TRUE,
                      zerolinecolor = "grey50", zerolinewidth = 1.5),
        yaxis = list(title = "-log<sub>10</sub>(FDR)", zeroline = TRUE,
                      zerolinecolor = "grey50", zerolinewidth = 1.5),
        zaxis = list(title = "Average Expression (AveExpr)"),
        camera = list(eye = list(x = 1.8, y = 1.2, z = 1.2))
      ),
      legend = list(title = list(text = "<b>Regulation</b>"))
    ) %>%
    config(displaylogo = FALSE, scrollZoom = TRUE, doubleClick = "reset")
  
  out_file <- file.path(repo_root, "results", "v3", "figures",
                         paste0("Volcano_", pw_id, "_3D.html"))
  saveWidget(fig, file = out_file, selfcontained = FALSE,
             title = paste0("3D Volcano \u2014 ", pw_labels[pw_id]))
  message(sprintf("  Saved: %s", out_file))
}

# ═══════════════════════════════════════════════════════════════════
# PART 2: 3D PPI Network — todos os 106 genes do metabolismo central
# ═══════════════════════════════════════════════════════════════════
message("\n=== 3D PPI Network ===")

# Carregar membership e DEGs
pw_membership <- rio::import(file.path(repo_root, "results", "tables", "pathway_gene_membership.csv"))
all_genes <- pw_membership$gene[pw_membership$in_matrix]

# Carregar DEG paired de cada via
deg_all <- do.call(rbind, lapply(names(pw_labels), function(pw_id) {
  f <- file.path(repo_root, "results", "v3", "tables", paste0("DEG_", pw_id, ".csv"))
  if (file.exists(f)) {
    d <- rio::import(f)
    d$pathway <- pw_id
    d
  }
}))
deg_all <- deg_all[!duplicated(deg_all$gene_id), ]

# Classificar regulação (paired, |logFC|>1, FDR<0.05)
deg_all$regulation <- "NS"
deg_all$regulation[deg_all$FDR_Paired < 0.05 & deg_all$logFC_Paired > 1] <- "Up"
deg_all$regulation[deg_all$FDR_Paired < 0.05 & deg_all$logFC_Paired < -1] <- "Down"

message(sprintf("  Genes to map: %d (%d Up, %d Down, %d NS)",
                nrow(deg_all),
                sum(deg_all$regulation == "Up", na.rm = TRUE),
                sum(deg_all$regulation == "Down", na.rm = TRUE),
                sum(deg_all$regulation == "NS", na.rm = TRUE)))

# ── STRING mapping ──
options(timeout = 5000000)
string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,
  score_threshold = 0,
  input_directory = ""
)

node_df <- data.frame(
  gene_symbol = deg_all$gene_id,
  logFC = deg_all$logFC_Paired,
  FDR = deg_all$FDR_Paired,
  regulation = deg_all$regulation,
  stringsAsFactors = FALSE
)

mapped <- string_db$map(node_df, "gene_symbol", removeUnmappedRows = TRUE)
mapped <- mapped[!is.na(mapped$STRING_id), ]
message(sprintf("  Mapped: %d/%d genes", nrow(mapped), nrow(node_df)))

if (nrow(mapped) < 2) {
  stop("STRING mapping returned fewer than 2 proteins. Cannot build network.")
}

# ── Get interactions ──
message("  Downloading STRING interactions...")
ppi_raw <- string_db$get_interactions(mapped$STRING_id)
score_cutoff <- 700

ppi <- ppi_raw %>%
  filter(from %in% mapped$STRING_id, to %in% mapped$STRING_id) %>%
  filter(combined_score >= score_cutoff)

message(sprintf("  Edges (combined_score >= %d): %d", score_cutoff, nrow(ppi)))

if (nrow(ppi) == 0) {
  stop("No PPI edges after high-confidence filtering.")
}

# ── Build graph ──
g <- graph_from_data_frame(
  d = ppi %>% transmute(from, to, weight = combined_score),
  directed = FALSE,
  vertices = mapped %>% transmute(
    name = STRING_id,
    gene_symbol = gene_symbol,
    regulation = regulation,
    logFC = logFC,
    FDR = FDR
  )
)

# Largest connected component
comp <- components(g)
giant_idx <- which.max(comp$csize)
g_cc <- induced_subgraph(g, vids = V(g)[comp$membership == giant_idx])
message(sprintf("  Largest component: %d nodes, %d edges", vcount(g_cc), ecount(g_cc)))

# ── 3D Layout ──
set.seed(42)
layout_3d <- layout_with_fr(g_cc, dim = 3, niter = 3000)

# ── Node attributes ──
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

# Escala de tamanho pelo degree
if (max(node_out$degree) > min(node_out$degree)) {
  node_out$size <- 4 + 10 * (node_out$degree - min(node_out$degree)) /
    (max(node_out$degree) - min(node_out$degree))
} else {
  node_out$size <- 8
}

# ── Build edge traces ──
el <- as_edgelist(g_cc)
edge_x <- c(); edge_y <- c(); edge_z <- c()

name_to_idx <- setNames(seq_len(nrow(node_out)), V(g_cc)$name)
for (i in seq_len(nrow(el))) {
  fi <- name_to_idx[el[i, 1]]
  ti <- name_to_idx[el[i, 2]]
  if (!is.na(fi) && !is.na(ti)) {
    edge_x <- c(edge_x, node_out$x[fi], node_out$x[ti], NA)
    edge_y <- c(edge_y, node_out$y[fi], node_out$y[ti], NA)
    edge_z <- c(edge_z, node_out$z[fi], node_out$z[ti], NA)
  }
}

# ── Build Plotly figure ──
fig_ppi <- plot_ly() %>%
  # Arestas
  add_trace(
    x = edge_x, y = edge_y, z = edge_z,
    type = "scatter3d", mode = "lines",
    line = list(color = "rgba(150,150,150,0.25)", width = 0.6),
    hoverinfo = "none", showlegend = FALSE
  )

# Nós por regulação
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
        size = ~size,
        color = col,
        opacity = op,
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
        "<b>PPI Network \u2014 Central Carbon Metabolism (3D)</b><br>",
        "<sup>STRING v11.5 | Comb. Score \u2265 ", score_cutoff,
        " | ", vcount(g_cc), " nodes | ", ecount(g_cc), " edges | ",
        "Layout: Fruchterman-Reingold 3D</sup>"
      ),
      font = list(size = 16)
    ),
    scene = list(
      xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE,
                    showticklabels = FALSE),
      yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE,
                    showticklabels = FALSE),
      zaxis = list(title = "", showgrid = FALSE, zeroline = FALSE,
                    showticklabels = FALSE),
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
    displaylogo = FALSE,
    scrollZoom = TRUE,
    doubleClick = "reset"
  )

out_ppi <- file.path(repo_root, "results", "v3", "figures", "PPI_network_3D.html")
saveWidget(fig_ppi, file = out_ppi, selfcontained = FALSE,
           title = "3D PPI Network \u2014 Central Carbon Metabolism in KIRP")
message(sprintf("  Saved: %s", out_ppi))

# Salvar dados do PPI
dir.create(file.path(repo_root, "results", "v3", "tables"), recursive = TRUE, showWarnings = FALSE)

# Tabela de centralidade 3D
rio::export(node_out, file.path(repo_root, "results", "v3", "tables", "ppi_3d_centrality.csv"))

# Resumo da rede
network_summary <- data.frame(
  metric = c("genes_input", "genes_mapped", "nodes_giant_component",
             "edges_giant_component", "score_cutoff", "n_up", "n_down", "n_ns"),
  value = c(nrow(node_df), nrow(mapped), vcount(g_cc), ecount(g_cc),
            score_cutoff,
            sum(node_out$regulation == "Up"),
            sum(node_out$regulation == "Down"),
            sum(node_out$regulation == "NS"))
)
rio::export(network_summary, file.path(repo_root, "results", "v3", "tables", "ppi_3d_summary.csv"))

message("\n=== 3D Visualizations Complete ===")
message(sprintf("Volcano 3D: results/v3/figures/Volcano_*_3D.html"))
message(sprintf("PPI 3D:    results/v3/figures/PPI_network_3D.html"))
