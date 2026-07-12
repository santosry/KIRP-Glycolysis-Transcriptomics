# 12_flowchart.R
# 🍌 Fluxograma analítico — Nano Banana Theme

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)

# ── 🍌 NANO BANANA PALETTE ──
banana_cream  <- "#FFF8DC"
banana_yellow <- "#FFE135"
banana_dark   <- "#D4A017"
banana_brown  <- "#8B6914"
banana_bg     <- "#FFFAF0"
banana_light  <- "#FFFDE7"
purple_node   <- "#8A2BE2"
blue_node     <- "#1F5BFF"
grey_node     <- "#E8D5A3"
node_text     <- "#3E2723"

# Build flowchart nodes
nodes <- tribble(
  ~id, ~label, ~x, ~y, ~phase,
  "N1", "UCSC Xena\nkidney.tsv", 1, 10, "data",
  "N2", "01 Data\nProvenance", 1, 8.5, "audit",
  "N3", "02 Prepare\n& QC", 1, 7, "qc",
  "N4", "03 Sample\nQC", 1, 5.5, "qc",
  
  "N5", "04 PCA\n& UMAP", 1, 4, "explore",
  
  "N6", "EIXO 1\nDEG Global", 4, 6, "global",
  "N7", "ORA KEGG\nUp | Down", 6.5, 7, "global",
  "N8", "ORA Reactome\nUp | Down", 6.5, 5, "global",
  
  "N9", "EIXO 2\nKEGG hsa00010", 1, 1.5, "targeted",
  "N10", "Heatmap\nhsa00010", 3, 0.5, "targeted",
  "N11", "Ranking\n|logFC|", 3, 2.5, "targeted",
  
  "N12", "STRING\nNetwork", 6.5, 3, "network",
  "N13", "Communities\nLeiden/Louvain", 6.5, 1.5, "network",
  "N14", "Centrality\nDegree/Between.", 6.5, 0, "network",
  
  "N15", "Integrative\nAnalysis", 9, 3, "integrate",
  "N16", "Sensitivity\nAnalysis", 9, 1.5, "integrate",
  "N17", "Manuscript\nResults", 9, 0, "output"
)

edges <- tribble(
  ~from, ~to,
  "N1", "N2", "N2", "N3", "N3", "N4",
  "N4", "N5",
  "N5", "N6",
  "N6", "N7", "N6", "N8",
  "N5", "N9",
  "N9", "N10", "N9", "N11",
  "N6", "N12",
  "N12", "N13", "N12", "N14",
  "N13", "N15", "N14", "N15",
  "N7", "N15", "N8", "N15",
  "N10", "N15", "N11", "N15",
  "N15", "N16",
  "N16", "N17"
)

# Build line segments for edges
edge_lines <- lapply(1:(length(edges$from)), function(i) {
  f <- edges$from[i]; t <- edges$to[i]
  nf <- nodes[nodes$id == f, ]; nt <- nodes[nodes$id == t, ]
  tibble(x = nf$x, xend = nt$x, y = nf$y, yend = nt$y)
}) |> bind_rows()

phase_colors <- c(
  data      = "#D4A017",
  audit     = "#C49A1C",
  qc        = "#B8942A",
  explore   = "#AB8E35",
  global    = blue_node,
  targeted  = purple_node,
  network   = "#E8963A",
  integrate = "#D4A017",
  output    = banana_brown
)

p_flow <- ggplot() +
  # Background
  annotate("rect", xmin = -0.5, xmax = 10.5, ymin = -1.5, ymax = 11,
           fill = banana_bg, alpha = 0.5) +
  # Edges
  geom_curve(data = edge_lines,
             aes(x = x, y = y, xend = xend, yend = yend),
             curvature = 0.1, color = "#C4A035", linewidth = 0.5, alpha = 0.6) +
  # Phase bands
  annotate("rect", xmin = -0.3, xmax = 2.5, ymin = 6.5, ymax = 10.5,
           fill = banana_cream, alpha = 0.4) +
  annotate("text", x = 1.1, y = 10.2, label = "DATA & QC", hjust = 0,
           color = banana_brown, fontface = "bold", size = 3.5) +
  annotate("rect", xmin = 3.5, xmax = 8, ymin = 4.5, ymax = 7.8,
           fill = "#E8E0F0", alpha = 0.3) +
  annotate("text", x = 3.7, y = 7.6, label = "EIXO 1: TRANSCRIPTOMA GLOBAL", hjust = 0,
           color = blue_node, fontface = "bold", size = 3.5) +
  annotate("rect", xmin = -0.3, xmax = 3.8, ymin = -0.5, ymax = 3.2,
           fill = "#F0E0F0", alpha = 0.3) +
  annotate("text", x = -0.1, y = 3.0, label = "EIXO 2: VIA hsa00010", hjust = 0,
           color = purple_node, fontface = "bold", size = 3.5) +
  annotate("rect", xmin = 5.8, xmax = 8, ymin = -0.8, ymax = 3.8,
           fill = "#FFF0E0", alpha = 0.3) +
  annotate("text", x = 6, y = 3.6, label = "REDE", hjust = 0,
           color = "#E8963A", fontface = "bold", size = 3.5) +
  annotate("rect", xmin = 8.3, xmax = 10.2, ymin = -0.8, ymax = 3.8,
           fill = banana_cream, alpha = 0.4) +
  annotate("text", x = 8.5, y = 3.6, label = "INTEGRAÇÃO", hjust = 0,
           color = banana_dark, fontface = "bold", size = 3.5) +
  # Title
  annotate("text", x = 5, y = 10.8,
           label = "🍌 PIPELINE TRANSCRIPTÔMICO KIRP — NANO BANANA 🍌",
           size = 6, fontface = "bold", color = banana_brown) +
  # Nodes
  geom_label(data = nodes, aes(x = x, y = y, label = label, fill = phase),
             color = node_text, size = 2.8, fontface = "bold",
             label.padding = unit(0.3, "lines"), label.r = unit(0.3, "lines"),
             alpha = 0.9) +
  scale_fill_manual(values = phase_colors, guide = "none") +
  xlim(-0.5, 10.5) + ylim(-1.5, 11) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = "grey80", linewidth = 1),
        plot.margin = margin(10, 10, 10, 10))

ggsave(file.path(repo_root, "results", "figures", "Flowchart_pipeline.png"),
       p_flow, width = 14, height = 11, dpi = 300)

message("🍌 Flowchart saved! Nano banana complete.")
