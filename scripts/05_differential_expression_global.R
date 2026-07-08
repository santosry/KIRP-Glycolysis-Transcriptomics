# 05_differential_expression_global.R
# EIXO 1: DEG GLOBAL com diagnóstico limma (plotSA, trend vs no-trend)

suppressPackageStartupMessages({
  library(limma); library(dplyr); library(tibble); library(rio); library(ggplot2)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
meta_file <- file.path(repo_root, "data", "processed", "metadata.rds")
expr_file <- file.path(repo_root, "data", "processed", "expression_matrix.rds")
if (!file.exists(meta_file) || !file.exists(expr_file)) stop("Run 02_prepare_data.R first.")

dir.create(file.path(repo_root, "results", "differential_expression"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

meta <- readRDS(meta_file); expr <- readRDS(expr_file)
E <- t(expr); storage.mode(E) <- "numeric"
gene_vars <- apply(E, 1, var, na.rm = TRUE)
E <- E[gene_vars > 0, ]
n_universe <- nrow(E)
message("Universe: ", n_universe, " genes")

# ── Design ──
design <- model.matrix(~ 0 + condition, data = meta)
colnames(design) <- gsub("condition", "", colnames(design))
stopifnot(ncol(E) == nrow(design))

fit <- lmFit(E, design)
contrast_matrix <- makeContrasts(KIRP_vs_Normal = KIRP - Normal_GTEx, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)

# ── DIAGNOSTIC: mean-variance trend ──
# Fit with and without trend, compare
fit2_trend <- eBayes(fit2, robust = TRUE, trend = TRUE)
fit2_notrend <- eBayes(fit2, robust = TRUE, trend = FALSE)

# plotSA
png(file.path(repo_root, "results", "figures", "limma_plotSA_trend.png"), width = 8, height = 6, units = "in", res = 200)
plotSA(fit2_trend, main = "Mean-Variance Trend (eBayes robust=TRUE, trend=TRUE)")
dev.off()

png(file.path(repo_root, "results", "figures", "limma_plotSA_notrend.png"), width = 8, height = 6, units = "in", res = 200)
plotSA(fit2_notrend, main = "Mean-Variance Trend (eBayes robust=TRUE, trend=FALSE)")
dev.off()

# Compare DEG counts between trend and no-trend
deg_trend <- topTable(fit2_trend, number = Inf, adjust.method = "BH", sort.by = "P")
deg_notrend <- topTable(fit2_notrend, number = Inf, adjust.method = "BH", sort.by = "P")

classify <- function(df, fdr = 0.05, lfc = 1) {
  df$regulation <- "NS"
  df$regulation[df$adj.P.Val < fdr & df$logFC > lfc] <- "Up"
  df$regulation[df$adj.P.Val < fdr & df$logFC < -lfc] <- "Down"
  df
}

dg1 <- classify(deg_trend); dg2 <- classify(deg_notrend)
message(sprintf("DEGs (trend=TRUE):  %d Up | %d Down", sum(dg1$regulation=="Up"), sum(dg1$regulation=="Down")))
message(sprintf("DEGs (trend=FALSE): %d Up | %d Down", sum(dg2$regulation=="Up"), sum(dg2$regulation=="Down")))

# Correlation between models
common_genes <- intersect(rownames(dg1), rownames(dg2))
logFC_cor <- cor(dg1[common_genes, "logFC"], dg2[common_genes, "logFC"])
dir_agree <- mean(sign(dg1[common_genes, "logFC"]) == sign(dg2[common_genes, "logFC"]))
message(sprintf("logFC correlation (trend vs no-trend): %.4f | Directional agreement: %.1f%%", logFC_cor, 100*dir_agree))

# ── USE trend=TRUE as primary (accounts for mean-variance relationship) ──
deg_all <- deg_trend |> as.data.frame() |> rownames_to_column("gene_id") |> as_tibble()

# ── Sensitivity: 3 thresholds ──
thresholds <- list(
  "FDR005_LFC0"  = list(fdr = 0.05, lfc = 0),
  "FDR005_LFC05" = list(fdr = 0.05, lfc = 0.5),
  "FDR005_LFC1"  = list(fdr = 0.05, lfc = 1.0)
)

deg_summaries <- list()
for (nm in names(thresholds)) {
  th <- thresholds[[nm]]
  dg <- deg_all
  dg$regulation <- "NS"
  dg$regulation[dg$adj.P.Val < th$fdr & dg$logFC > th$lfc] <- "Up"
  dg$regulation[dg$adj.P.Val < th$fdr & dg$logFC < -th$lfc] <- "Down"
  
  n_u <- sum(dg$regulation == "Up"); n_d <- sum(dg$regulation == "Down")
  deg_summaries[[nm]] <- tibble(threshold=nm, fdr=th$fdr, lfc=th$lfc, n_up=n_u, n_down=n_d, n_total=n_u+n_d)
  
  if (nm == "FDR005_LFC1") deg_primary <- dg
}
deg_sens <- bind_rows(deg_summaries)
rio::export(deg_sens, file.path(repo_root, "results", "tables", "deg_sensitivity.csv"))

# ── Save primary ──
rio::export(deg_primary, file.path(repo_root, "results", "differential_expression", "DEG_global.csv"))
rio::export(tibble(gene_id=rownames(E)), file.path(repo_root, "results", "differential_expression", "gene_universe.csv"))

# ── Volcano plot ──
deg_primary$color <- deg_primary$regulation
p_volc <- ggplot(deg_primary, aes(logFC, -log10(adj.P.Val), color=color)) +
  geom_point(size=0.6, alpha=0.5) +
  scale_color_manual(values=c(Up="#1F5BFF", Down="#8A2BE2", NS="grey70")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", linewidth=0.4) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth=0.4) +
  labs(title="Volcano Plot — Global DE",
       subtitle=paste0(sum(deg_primary$regulation=="Up")," Up | ",sum(deg_primary$regulation=="Down")," Down"),
       x="log2 Fold Change", y="-log10(FDR)") + theme_classic(14)
ggsave(file.path(repo_root, "results", "figures", "Volcano_global.png"), p_volc, width=8, height=7, dpi=300)

# ── MA plot ──
p_ma <- ggplot(deg_primary, aes(AveExpr, logFC, color=color)) +
  geom_point(size=0.6, alpha=0.5) +
  scale_color_manual(values=c(Up="#1F5BFF", Down="#8A2BE2", NS="grey70")) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.4) +
  labs(title="MA Plot — Global DE", x="Average Expression", y="log2 Fold Change") + theme_classic(14)
ggsave(file.path(repo_root, "results", "figures", "MA_plot_global.png"), p_ma, width=8, height=6, dpi=300)

# ── Diagnostic summary ──
diag_summary <- tibble(
  metric = c("n_genes_tested", "n_DEG_trend_Up", "n_DEG_trend_Down", "n_DEG_notrend_Up", "n_DEG_notrend_Down",
             "logFC_correlation", "directional_agreement"),
  value = c(n_universe, sum(dg1$regulation=="Up"), sum(dg1$regulation=="Down"),
            sum(dg2$regulation=="Up"), sum(dg2$regulation=="Down"),
            round(logFC_cor, 4), round(dir_agree, 4))
)
rio::export(diag_summary, file.path(repo_root, "results", "tables", "limma_diagnostics.csv"))

message("\n✓ Global DE complete: ", sum(deg_primary$regulation=="Up"), " Up | ", sum(deg_primary$regulation=="Down"), " Down")
