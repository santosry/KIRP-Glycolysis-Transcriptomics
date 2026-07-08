# 05_differential_expression_global.R
# EIXO 1: ANÁLISE TRANSCRIPTÔMICA GLOBAL COM SENSIBILIDADE

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tibble)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
meta_file <- file.path(repo_root, "data", "processed", "metadata.rds")
expr_file <- file.path(repo_root, "data", "processed", "expression_matrix.rds")

if (!file.exists(meta_file) || !file.exists(expr_file)) {
  stop("Run 02_prepare_data.R first.")
}

dir.create(file.path(repo_root, "results", "differential_expression"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)

meta <- readRDS(meta_file)
expr <- readRDS(expr_file)

# Transpose: samples × genes → genes × samples
E <- t(expr)
storage.mode(E) <- "numeric"

# Remove constant genes
gene_vars <- apply(E, 1, var, na.rm = TRUE)
E <- E[gene_vars > 0, ]
n_universe <- nrow(E)
message("Universe of genes tested: ", n_universe)

# ── Model ──
# Input is log2(norm_count + 1) from UCSC Xena → limma with eBayes() is appropriate.
# For raw counts, voom() would be required. The data is pre-normalized.
# eBayes borrows information across genes, stabilizing variance estimates
# especially important given the sample imbalance (KIRP >> Normal).

design <- model.matrix(~ 0 + condition, data = meta)
colnames(design) <- levels(meta$condition)
stopifnot(ncol(E) == nrow(design))

fit <- lmFit(E, design)
contrast_matrix <- makeContrasts(KIRP_vs_Normal = KIRP - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, robust = TRUE)

# ── Full results table ──
deg_all <- topTable(fit2, coef = "KIRP_vs_Normal", number = Inf,
                     adjust.method = "BH", sort.by = "P") |>
  rownames_to_column(var = "gene_id") |>
  as_tibble()

# ── Classify DEGs at multiple thresholds ──
thresholds <- list(
  "FDR005_logFC0"   = list(fdr = 0.05, lfc = 0),
  "FDR005_logFC05"  = list(fdr = 0.05, lfc = 0.5),
  "FDR005_logFC1"   = list(fdr = 0.05, lfc = 1.0)
)

deg_tables <- list()
deg_summaries <- list()

for (thresh_name in names(thresholds)) {
  th <- thresholds[[thresh_name]]
  
  deg_copy <- deg_all
  deg_copy$regulation <- "NS"
  deg_copy$regulation[deg_copy$adj.P.Val < th$fdr & deg_copy$logFC > th$lfc] <- "Up"
  deg_copy$regulation[deg_copy$adj.P.Val < th$fdr & deg_copy$logFC < -th$lfc] <- "Down"
  
  deg_tables[[thresh_name]] <- deg_copy
  
  n_up <- sum(deg_copy$regulation == "Up")
  n_down <- sum(deg_copy$regulation == "Down")
  
  deg_summaries[[thresh_name]] <- tibble(
    threshold = thresh_name,
    fdr_cutoff = th$fdr,
    lfc_cutoff = th$lfc,
    n_up = n_up,
    n_down = n_down,
    n_total_deg = n_up + n_down
  )
  
  message(sprintf("  %s: %d Up | %d Down | %d total DEGs", thresh_name, n_up, n_down, n_up + n_down))
}

# ── Save primary results (FDR < 0.05, |logFC| > 1) ──
deg_primary <- deg_tables[["FDR005_logFC1"]]
rio::export(deg_primary, file.path(repo_root, "results", "differential_expression", "DEG_global.csv"))

# ── Save sensitivity comparison ──
deg_sensitivity <- bind_rows(deg_summaries)
rio::export(deg_sensitivity, file.path(repo_root, "results", "tables", "deg_sensitivity.csv"))

# ── Export universe gene list ──
universe_genes <- tibble(gene_id = rownames(E))
rio::export(universe_genes, file.path(repo_root, "results", "differential_expression", "gene_universe.csv"))

# ── MA plot ──
library(ggplot2)
deg_primary$color <- "NS"
deg_primary$color[deg_primary$regulation == "Up"] <- "Up"
deg_primary$color[deg_primary$regulation == "Down"] <- "Down"

p_ma <- ggplot(deg_primary, aes(x = AveExpr, y = logFC, color = color)) +
  geom_point(size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c(Up = "#1F5BFF", Down = "#8A2BE2", NS = "grey70")) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  labs(title = "MA plot — Expressão diferencial global",
       x = "Expressão média (AveExpr)", y = "log2 Fold Change", color = "Regulação") +
  theme_classic(base_size = 14)

ggsave(file.path(repo_root, "results", "figures", "MA_plot_global.png"), p_ma, width = 8, height = 6, dpi = 300)

# ── Volcano plot global with DEG counts ──
p_volcano <- ggplot(deg_primary, aes(x = logFC, y = -log10(adj.P.Val), color = color)) +
  geom_point(size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c(Up = "#1F5BFF", Down = "#8A2BE2", NS = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4, color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, color = "grey40") +
  labs(title = "Volcano plot — Expressão diferencial global",
       subtitle = paste0(sum(deg_primary$regulation == "Up"), " Up | ",
                        sum(deg_primary$regulation == "Down"), " Down | ",
                        sum(deg_primary$regulation == "NS"), " NS"),
       x = "log2 Fold Change", y = "-log10(FDR)", color = "Regulação") +
  theme_classic(base_size = 14)

ggsave(file.path(repo_root, "results", "figures", "Volcano_global.png"), p_volcano, width = 8, height = 7, dpi = 300)

message("\n✓ Global DE complete. Universe: ", n_universe, " genes tested.")
