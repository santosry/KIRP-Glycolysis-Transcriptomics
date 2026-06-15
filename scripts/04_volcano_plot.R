suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
deg_file <- file.path(repo_root, "results", "differential_expression", "DEG_KIRP_vs_Normal.csv")

if (!file.exists(deg_file)) {
  stop("Run scripts/03_differential_expression.R before generating the volcano plot.")
}

dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)

deg <- rio::import(deg_file)
lfc_cutoff <- 1
fdr_cutoff <- 0.05

deg$volcano_class <- "NS"
deg$volcano_class[deg$adj.P.Val < fdr_cutoff & deg$logFC > lfc_cutoff] <- "Up"
deg$volcano_class[deg$adj.P.Val < fdr_cutoff & deg$logFC < -lfc_cutoff] <- "Down"

p <- ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val), color = volcano_class)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c(Up = "#1F5BFF", Down = "#8A2BE2", NS = "grey70")) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", linewidth = 0.4) +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", linewidth = 0.4) +
  geom_text_repel(
    data = subset(deg, volcano_class %in% c("Up", "Down")),
    aes(label = gene_symbol),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +
  labs(
    title = "Differential expression in glycolysis/gluconeogenesis genes",
    x = "log2 fold change",
    y = "-log10(FDR)",
    color = "Regulation"
  ) +
  theme_classic(base_size = 14)

ggsave(
  filename = file.path(repo_root, "results", "figures", "Volcano.png"),
  plot = p,
  width = 7,
  height = 6,
  dpi = 300
)

message("Saved results/figures/Volcano.png")
