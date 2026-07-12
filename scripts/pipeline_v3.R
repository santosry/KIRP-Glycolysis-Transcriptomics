# pipeline_v3.R — VERSÃO FINAL (AJUSTADA PARA kidney.tsv COM 110 GENES)
# Matriz contém apenas os genes das 3 vias KEGG
# Análise pareada primária, extração pós-hoc por via, concordância

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tibble)
  library(rio)
  library(ggplot2)
  library(ggrepel)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
setwd(repo_root)

dir.create("results/v3", recursive = TRUE, showWarnings = FALSE)
dir.create("results/v3/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/v3/tables", recursive = TRUE, showWarnings = FALSE)

# ═════════════════════════════════════════
# STEP 1: Load data
# ═════════════════════════════════════════
message("=== STEP 1: Loading kidney.tsv ===")
raw <- rio::import("data/raw/kidney.tsv")
n_genes_raw <- ncol(raw) - 8  # 8 metadata columns
n_samples <- nrow(raw)
message(sprintf("Raw: %d genes x %d samples", n_genes_raw, n_samples))

# Separate metadata and expression
meta_cols <- c("sample", "samples", "TCGA_GTEX_main_category", "_sample_type", 
               "_study", "_primary_site", "OS", "OS.time")
gene_symbols <- setdiff(colnames(raw), meta_cols)

expr_data <- as.matrix(raw[, gene_symbols])
rownames(expr_data) <- raw$sample
storage.mode(expr_data) <- "numeric"

# Transpose: genes x samples
E_all <- t(expr_data)
message(sprintf("Expression matrix: %d genes x %d samples", nrow(E_all), ncol(E_all)))

# Remove low-expression genes (expression <= 1 in >90% of samples)
n_low <- sum(rowMeans(E_all > 1) <= 0.1)
E_filt <- E_all[rowMeans(E_all > 1) > 0.1, ]
n_genes_final <- nrow(E_filt)
genes_removed <- setdiff(rownames(E_all), rownames(E_filt))
message(sprintf("Filtering: %d -> remove %d low expression = %d final", 
                nrow(E_all), n_low, n_genes_final))
if (length(genes_removed) > 0) {
  message(sprintf("  Removed: %s", paste(genes_removed, collapse=", ")))
}

# ═════════════════════════════════════════
# STEP 2: Metadata
# ═════════════════════════════════════════
message("\n=== STEP 2: Metadata ===")
samples <- colnames(E_filt)
meta <- tibble(
  sample = raw$sample,
  study = ifelse(grepl("^GTEX-", raw$sample), "GTEX", "TCGA"),
  sample_type = case_when(
    grepl("-01", raw$sample, fixed = TRUE) ~ "Primary_Tumor",
    grepl("-11", raw$sample, fixed = TRUE) ~ "Solid_Tissue_Normal",
    grepl("^GTEX-", raw$sample) ~ "Normal_Tissue",
    TRUE ~ "Other"
  ),
  condition = case_when(
    sample_type == "Primary_Tumor" ~ "KIRP",
    study == "GTEX" ~ "Normal_GTEx",
    sample_type == "Solid_Tissue_Normal" ~ "TCGA_Normal",
    TRUE ~ "Other"
  ),
  participant = sub("^(TCGA-..-....).*", "\\1", raw$sample)
)

kirp_participants <- unique(meta$participant[meta$condition == "KIRP"])
meta$paired_normal <- meta$condition == "TCGA_Normal" & meta$participant %in% kirp_participants

message("Sample counts:")
print(table(meta$condition))
message(sprintf("TCGA Normals: %d (KIRP-adjacent: %d, Other: %d)", 
                sum(meta$condition == "TCGA_Normal"),
                sum(meta$paired_normal),
                sum(meta$condition == "TCGA_Normal" & !meta$paired_normal)))

# Restrict metadata to samples in expression matrix
meta <- meta[meta$sample %in% colnames(E_filt), ]

# ═════════════════════════════════════════
# STEP 3: QC — PCA
# ═════════════════════════════════════════
message("\n=== QC: PCA ===")
pca <- prcomp(t(E_filt), scale. = TRUE, center = TRUE)
pca_scores <- as.data.frame(pca$x)
pca_scores$condition <- meta$condition[match(rownames(pca_scores), meta$sample)]
pca_scores$study <- meta$study[match(rownames(pca_scores), meta$sample)]
var_pc1 <- round(summary(pca)$importance[2,1] * 100, 1)
var_pc2 <- round(summary(pca)$importance[2,2] * 100, 1)

png("results/v3/figures/PCA_transcriptome.png", width = 1000, height = 800, res = 130)
ggplot(pca_scores, aes(PC1, PC2, color = condition, shape = study)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(title = "PCA — 110 Central Carbon Metabolism Genes",
       subtitle = sprintf("PC1: %.1f%% | PC2: %.1f%% | Colored by condition, shaped by study", var_pc1, var_pc2)) +
  theme_minimal(14)
dev.off()

# QC density
png("results/v3/figures/QC_density.png", width = 900, height = 600, res = 120)
kirp_cols <- meta$sample[meta$condition == "KIRP"]
tcga_cols <- meta$sample[meta$condition == "TCGA_Normal"]
gtex_cols <- meta$sample[meta$condition == "Normal_GTEx"]
plot(density(as.vector(E_filt[, kirp_cols[1:min(50, length(kirp_cols))]])), 
     col = "#8A2BE2", lwd = 2, main = "Expression Density by Condition", 
     xlab = "log2(expression)", ylim = c(0, 0.25))
lines(density(as.vector(E_filt[, tcga_cols[1:min(50, length(tcga_cols))]])), col = "#1F5BFF", lwd = 2)
lines(density(as.vector(E_filt[, gtex_cols])), col = "#FFB347", lwd = 2)
legend("topright", c("KIRP", "TCGA Normal", "GTEx"), col = c("#8A2BE2", "#1F5BFF", "#FFB347"), lwd = 2)
dev.off()

# ═════════════════════════════════════════
# STEP 4: PRIMARY — Paired analysis (32 pairs)
# ═════════════════════════════════════════
message("\n=== STEP 4: PRIMARY — Paired (32 pairs) ===")
paired_normals <- meta$sample[meta$paired_normal]
paired_tumors <- meta$sample[meta$condition == "KIRP" & 
                              meta$participant %in% meta$participant[meta$paired_normal]]
paired_meta <- meta[meta$sample %in% c(paired_normals, paired_tumors), ]
paired_meta$condition <- factor(paired_meta$condition, levels = c("TCGA_Normal", "KIRP"))
paired_meta$patient <- factor(paired_meta$participant)

E_paired <- E_filt[, paired_meta$sample, drop = FALSE]
design_paired <- model.matrix(~ patient + condition, data = paired_meta)
fit_paired <- lmFit(E_paired, design_paired)
fit_paired <- eBayes(fit_paired, robust = TRUE, trend = TRUE)

coef_idx <- grep("conditionKIRP", colnames(design_paired))
deg_paired <- topTable(fit_paired, coef = coef_idx, number = Inf, adjust.method = "BH", confint = TRUE)
deg_paired$gene_id <- rownames(deg_paired)
deg_paired$regulation <- ifelse(deg_paired$adj.P.Val < 0.05 & abs(deg_paired$logFC) > 1,
                                 ifelse(deg_paired$logFC > 0, "Up", "Down"), "NS")
n_up <- sum(deg_paired$regulation == "Up")
n_down <- sum(deg_paired$regulation == "Down")
message(sprintf("Paired DEGs: %d Up | %d Down | %d NS (of %d genes)", 
                n_up, n_down, sum(deg_paired$regulation=="NS"), nrow(deg_paired)))

png("results/v3/figures/plotSA_paired.png", width = 800, height = 600, res = 120)
plotSA(fit_paired, main = sprintf("Mean-Variance — Paired Analysis (%d genes)", nrow(E_filt)))
dev.off()

# ═════════════════════════════════════════
# STEP 5: SECONDARY — KIRP vs KIRP-adjacent
# ═════════════════════════════════════════
message("\n=== STEP 5: SECONDARY — KIRP vs Adjacent ===")
tcga_idx <- meta$condition == "KIRP" | meta$paired_normal
tcga_meta <- meta[tcga_idx, ]
tcga_meta$condition <- factor(tcga_meta$condition, levels = c("TCGA_Normal", "KIRP"))
E_tcga <- E_filt[, tcga_meta$sample, drop = FALSE]

design_tcga <- model.matrix(~ condition, data = tcga_meta)
fit_tcga <- lmFit(E_tcga, design_tcga)
fit_tcga <- eBayes(fit_tcga, robust = TRUE, trend = TRUE)

deg_tcga <- topTable(fit_tcga, coef = "conditionKIRP", number = Inf, adjust.method = "BH", confint = TRUE)
deg_tcga$gene_id <- rownames(deg_tcga)
message(sprintf("KIRP vs adjacent: %d Up | %d Down", 
                sum(deg_tcga$adj.P.Val < 0.05 & deg_tcga$logFC > 1),
                sum(deg_tcga$adj.P.Val < 0.05 & deg_tcga$logFC < -1)))

# ═════════════════════════════════════════
# STEP 6: EXPLORATORY — KIRP vs GTEx
# ═════════════════════════════════════════
message("\n=== STEP 6: EXPLORATORY — KIRP vs GTEx ===")
gtx_idx <- meta$condition %in% c("KIRP", "Normal_GTEx")
gtx_meta <- meta[gtx_idx, ]
gtx_meta$condition <- factor(gtx_meta$condition, levels = c("Normal_GTEx", "KIRP"))
E_gtx <- E_filt[, gtx_meta$sample, drop = FALSE]

design_gtx <- model.matrix(~ condition, data = gtx_meta)
fit_gtx <- lmFit(E_gtx, design_gtx)
fit_gtx <- eBayes(fit_gtx, robust = TRUE, trend = TRUE)

deg_gtx <- topTable(fit_gtx, coef = "conditionKIRP", number = Inf, adjust.method = "BH", confint = TRUE)
deg_gtx$gene_id <- rownames(deg_gtx)
message(sprintf("KIRP vs GTEx: %d Up | %d Down (WARNING: confounded)",
                sum(deg_gtx$adj.P.Val < 0.05 & deg_gtx$logFC > 1),
                sum(deg_gtx$adj.P.Val < 0.05 & deg_gtx$logFC < -1)))

# ═════════════════════════════════════════
# STEP 7: Paired plots for key genes
# ═════════════════════════════════════════
message("\n=== STEP 7: Paired plots ===")
key_genes <- intersect(c("ALDOB", "HK2", "PCK1", "G6PD", "TKT", "ADH1B", "FBP1"), rownames(E_filt))
for (g in key_genes) {
  pd <- paired_meta
  pd$expression <- E_filt[g, pd$sample]
  pd$pair_id <- pd$participant
  
  png(sprintf("results/v3/figures/Paired_%s.png", g), width = 700, height = 500, res = 110)
  p <- ggplot(pd, aes(condition, expression, group = pair_id)) +
    geom_line(alpha = 0.4, color = "grey50") +
    geom_point(aes(color = condition), size = 3, alpha = 0.8) +
    scale_color_manual(values = c(TCGA_Normal = "#1F5BFF", KIRP = "#8A2BE2")) +
    labs(title = g, y = "log2(expression)", x = "") +
    theme_minimal(14) + theme(legend.position = "none")
  print(p)
  dev.off()
}

# ═════════════════════════════════════════
# STEP 8: Gene set testing (camera)
# ═════════════════════════════════════════
message("\n=== STEP 8: Gene set testing (camera) ===")
pw_membership <- rio::import("results/tables/pathway_gene_membership.csv")

pathways <- list(
  hsa00010_Glycolysis = intersect(pw_membership$gene[pw_membership$in_hsa00010], rownames(E_filt)),
  hsa00030_PPP        = intersect(pw_membership$gene[pw_membership$in_hsa00030], rownames(E_filt)),
  hsa00020_TCA        = intersect(pw_membership$gene[pw_membership$in_hsa00020], rownames(E_filt))
)

camera_res <- camera(E_paired, pathways, design_paired, coef = coef_idx)
camera_res$pathway <- rownames(camera_res)
message("Camera gene set test (paired):")
print(camera_res)
rio::export(camera_res, "results/v3/tables/camera_gene_sets.csv")

# ═════════════════════════════════════════
# STEP 9: Pathway extraction & DEG tables
# ═════════════════════════════════════════
message("\n=== STEP 9: Pathway extraction ===")
pw_names <- list(
  hsa00010 = "Glycolysis / Gluconeogenesis",
  hsa00030 = "Pentose Phosphate Pathway",
  hsa00020 = "Citrate Cycle (TCA)"
)

pw_results <- list()
for (pw_id in names(pw_names)) {
  pw_genes <- intersect(pw_membership$gene[pw_membership[[paste0("in_", pw_id)]]], rownames(E_filt))
  
  pw_paired <- deg_paired[deg_paired$gene_id %in% pw_genes, ]
  pw_tcga <- deg_tcga[deg_tcga$gene_id %in% pw_genes, ]
  pw_gtx <- deg_gtx[deg_gtx$gene_id %in% pw_genes, ]
  
  pw_out <- data.frame(
    gene_id = sort(pw_genes),
    stringsAsFactors = FALSE
  )
  pw_out$logFC_Paired <- pw_paired$logFC[match(pw_out$gene_id, pw_paired$gene_id)]
  pw_out$CI.L_Paired   <- pw_paired$CI.L[match(pw_out$gene_id, pw_paired$gene_id)]
  pw_out$CI.R_Paired   <- pw_paired$CI.R[match(pw_out$gene_id, pw_paired$gene_id)]
  pw_out$FDR_Paired    <- pw_paired$adj.P.Val[match(pw_out$gene_id, pw_paired$gene_id)]
  pw_out$AveExpr       <- pw_paired$AveExpr[match(pw_out$gene_id, pw_paired$gene_id)]
  pw_out$logFC_TCGA_adj <- pw_tcga$logFC[match(pw_out$gene_id, pw_tcga$gene_id)]
  pw_out$FDR_TCGA_adj  <- pw_tcga$adj.P.Val[match(pw_out$gene_id, pw_tcga$gene_id)]
  pw_out$logFC_GTEx    <- pw_gtx$logFC[match(pw_out$gene_id, pw_gtx$gene_id)]
  pw_out$FDR_GTEx      <- pw_gtx$adj.P.Val[match(pw_out$gene_id, pw_gtx$gene_id)]
  pw_out$pathway <- pw_id
  
  rio::export(pw_out, sprintf("results/v3/tables/DEG_%s.csv", pw_id))
  pw_results[[pw_id]] <- pw_out
  
  n_up <- sum(pw_out$FDR_Paired < 0.05 & pw_out$logFC_Paired > 1, na.rm = TRUE)
  n_down <- sum(pw_out$FDR_Paired < 0.05 & pw_out$logFC_Paired < -1, na.rm = TRUE)
  message(sprintf("  %s: %d genes, %d Up, %d Down", pw_id, nrow(pw_out), n_up, n_down))
}

# ═════════════════════════════════════════
# STEP 10: Supplementary Table S1
# ═════════════════════════════════════════
message("\n=== STEP 10: Supplementary table ===")
all_pw_genes <- unique(unlist(lapply(pw_results, `[[`, "gene_id")))
supp <- data.frame(gene_id = all_pw_genes, stringsAsFactors = FALSE)
supp$pathway <- ""
for (pw_id in names(pw_names)) {
  pw_genes <- pw_results[[pw_id]]$gene_id
  for (g in pw_genes) {
    if (supp$pathway[supp$gene_id == g] == "") {
      supp$pathway[supp$gene_id == g] <- pw_id
    } else {
      supp$pathway[supp$gene_id == g] <- paste(supp$pathway[supp$gene_id == g], pw_id, sep = ";")
    }
  }
}

supp$logFC_Paired <- deg_paired$logFC[match(supp$gene_id, deg_paired$gene_id)]
supp$CI.L_Paired   <- deg_paired$CI.L[match(supp$gene_id, deg_paired$gene_id)]
supp$CI.R_Paired   <- deg_paired$CI.R[match(supp$gene_id, deg_paired$gene_id)]
supp$P.Value_Paired <- deg_paired$P.Value[match(supp$gene_id, deg_paired$gene_id)]
supp$FDR_Paired    <- deg_paired$adj.P.Val[match(supp$gene_id, deg_paired$gene_id)]
supp$AveExpr       <- deg_paired$AveExpr[match(supp$gene_id, deg_paired$gene_id)]
supp$logFC_TCGA_adj <- deg_tcga$logFC[match(supp$gene_id, deg_tcga$gene_id)]
supp$FDR_TCGA_adj  <- deg_tcga$adj.P.Val[match(supp$gene_id, deg_tcga$gene_id)]
supp$logFC_GTEx    <- deg_gtx$logFC[match(supp$gene_id, deg_gtx$gene_id)]
supp$FDR_GTEx      <- deg_gtx$adj.P.Val[match(supp$gene_id, deg_gtx$gene_id)]

rio::export(supp, "results/v3/tables/Supplementary_Table_S1.csv")

n_deg <- sum(supp$FDR_Paired < 0.05 & abs(supp$logFC_Paired) > 1, na.rm = TRUE)
n_up_s <- sum(supp$FDR_Paired < 0.05 & supp$logFC_Paired > 1, na.rm = TRUE)
n_down_s <- sum(supp$FDR_Paired < 0.05 & supp$logFC_Paired < -1, na.rm = TRUE)
n_shared <- sum(supp$FDR_Paired < 0.05 & abs(supp$logFC_Paired) > 1 & grepl(";", supp$pathway), na.rm = TRUE)
message(sprintf("Unique DEGs: %d (%d Up, %d Down, %d shared)", n_deg, n_up_s, n_down_s, n_shared))

# ═════════════════════════════════════════
# STEP 11: Concordance
# ═════════════════════════════════════════
message("\n=== STEP 11: Concordance ===")

calc_concordance <- function(x, y, label_x, label_y, gene_ids) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]; ids <- gene_ids[ok]
  
  mx <- mean(x); my <- mean(y)
  vx <- var(x); vy <- var(y)
  sxy <- cov(x, y)
  ccc <- 2 * sxy / (vx + vy + (mx - my)^2)
  r <- cor(x, y)
  mae <- mean(abs(x - y))
  bias <- mean(x - y)
  
  lm_fit <- lm(y ~ x)
  slope <- coef(lm_fit)[2]
  intercept <- coef(lm_fit)[1]
  
  substantial <- abs(x) > 0.5 | abs(y) > 0.5
  dir_agree <- mean(sign(x[substantial]) == sign(y[substantial]))
  
  discordant_idx <- which(sign(x) != sign(y))
  
  list(ccc = ccc, r = r, mae = mae, bias = bias, slope = slope, intercept = intercept,
       dir_agree = dir_agree, n_total = length(x), n_substantial = sum(substantial),
       n_discordant = length(discordant_idx),
       label_x = label_x, label_y = label_y, 
       discordant_genes = ids[discordant_idx])
}

# Paired vs TCGA-adjacent
conc1 <- calc_concordance(supp$logFC_Paired, supp$logFC_TCGA_adj,
                           "Paired", "KIRP vs Adjacent", supp$gene_id)
message(sprintf("\nPaired vs TCGA-adjacent:"))
message(sprintf("  CCC = %.3f | r = %.3f | MAE = %.3f | Bias = %.3f", conc1$ccc, conc1$r, conc1$mae, conc1$bias))
message(sprintf("  Directional (>|0.5|): %.1f%% (%d/%d) | Discordant: %d", 
                conc1$dir_agree*100, sum(abs(supp$logFC_Paired)>0.5 | abs(supp$logFC_TCGA_adj)>0.5, na.rm=TRUE),
                sum(abs(supp$logFC_Paired)>0.5 | abs(supp$logFC_TCGA_adj)>0.5, na.rm=TRUE), conc1$n_discordant))

# Paired vs GTEx
conc2 <- calc_concordance(supp$logFC_Paired, supp$logFC_GTEx,
                           "Paired", "KIRP vs GTEx", supp$gene_id)
message(sprintf("\nPaired vs GTEx:"))
message(sprintf("  CCC = %.3f | r = %.3f | MAE = %.3f | Bias = %.3f", conc2$ccc, conc2$r, conc2$mae, conc2$bias))
message(sprintf("  Directional (>|0.5|): %.1f%% | Discordant: %d", conc2$dir_agree*100, conc2$n_discordant))

# Discordant genes table
discordant_df <- data.frame(
  gene_id = conc2$discordant_genes,
  logFC_Paired = supp$logFC_Paired[match(conc2$discordant_genes, supp$gene_id)],
  logFC_GTEx = supp$logFC_GTEx[match(conc2$discordant_genes, supp$gene_id)],
  FDR_Paired = supp$FDR_Paired[match(conc2$discordant_genes, supp$gene_id)],
  FDR_GTEx = supp$FDR_GTEx[match(conc2$discordant_genes, supp$gene_id)]
)
discordant_df <- discordant_df[order(abs(discordant_df$logFC_Paired), decreasing = TRUE), ]
rio::export(discordant_df, "results/v3/tables/discordant_genes.csv")

# Bland-Altman
ba_data <- data.frame(
  mean_logFC = (supp$logFC_Paired + supp$logFC_GTEx) / 2,
  diff_logFC = supp$logFC_Paired - supp$logFC_GTEx
)
ba_data <- ba_data[complete.cases(ba_data), ]

png("results/v3/figures/BlandAltman_Paired_vs_GTEx.png", width = 900, height = 700, res = 120)
ggplot(ba_data, aes(mean_logFC, diff_logFC)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_hline(yintercept = 0, linetype = "solid") +
  geom_hline(yintercept = conc2$bias, linetype = "dashed", color = "#8A2BE2") +
  geom_hline(yintercept = conc2$bias + 1.96 * sd(ba_data$diff_logFC), linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = conc2$bias - 1.96 * sd(ba_data$diff_logFC), linetype = "dotted", color = "grey50") +
  labs(title = "Bland-Altman: Paired vs GTEx",
       subtitle = sprintf("Bias = %.2f | LoA = [%.2f, %.2f]", 
                          conc2$bias, conc2$bias - 1.96*sd(ba_data$diff_logFC),
                          conc2$bias + 1.96*sd(ba_data$diff_logFC)),
       x = "Mean log2FC", y = "Difference (Paired - GTEx)") +
  theme_minimal(14)
dev.off()

# ═════════════════════════════════════════
# STEP 12: Volcano plots
# ═════════════════════════════════════════
message("\n=== STEP 12: Volcano plots ===")
for (pw_id in names(pw_names)) {
  pw_out <- pw_results[[pw_id]]
  pw_out <- pw_out[order(abs(pw_out$logFC_Paired), decreasing = TRUE), ]
  pw_out$regulation <- ifelse(pw_out$FDR_Paired < 0.05 & abs(pw_out$logFC_Paired) > 1,
                               ifelse(pw_out$logFC_Paired > 0, "Up", "Down"), "NS")
  pw_out$label <- ifelse(pw_out$regulation != "NS" | abs(pw_out$logFC_Paired) > 1.5, pw_out$gene_id, "")
  
  n_up <- sum(pw_out$regulation == "Up", na.rm = TRUE)
  n_down <- sum(pw_out$regulation == "Down", na.rm = TRUE)
  
  p <- ggplot(pw_out, aes(logFC_Paired, -log10(FDR_Paired), color = regulation, label = label)) +
    geom_point(alpha = 0.7, size = 2.5) +
    scale_color_manual(values = c(Up = "#1F5BFF", Down = "#8A2BE2", NS = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) +
    geom_text_repel(size = 3, max.overlaps = 40, box.padding = 0.3) +
    labs(title = paste0(pw_names[[pw_id]], " (Paired, 32 pairs)"),
         subtitle = sprintf("%d genes | %d Up | %d Down", nrow(pw_out), n_up, n_down),
         x = "log2(Fold Change) KIRP vs Adjacent Normal", y = "-log10(FDR)") +
    theme_minimal(14) + theme(plot.background = element_rect(fill = "white", color = NA))
  
  ggsave(sprintf("results/v3/figures/Volcano_%s.png", pw_id), p, width = 10, height = 8, dpi = 300)
}

# ═════════════════════════════════════════
# STEP 13: Summary
# ═════════════════════════════════════════
message("\n=== PIPELINE V3 COMPLETE ===")
message(sprintf("Genes analyzed: %d (after low-expression filter)", n_genes_final))
message(sprintf("Primary: Paired analysis, 32 pairs, %d DEGs", n_deg))
message(sprintf("Camera: %s", paste(rownames(camera_res), collapse=", ")))

# Save session info
writeLines(capture.output(sessionInfo()), "results/v3/sessionInfo.txt")
message("Session info saved")
