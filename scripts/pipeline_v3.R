# pipeline_v3.R — VERSÃO FINAL
# Transcriptoma completo, análise pareada primária,
# limma único, extração pós-hoc, concordância com métricas apropriadas

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tibble)
  library(rio)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(org.Hs.eg.db)
  library(clusterProfiler)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
setwd(repo_root)

dir.create("results/v3", recursive = TRUE, showWarnings = FALSE)
dir.create("results/v3/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/v3/tables", recursive = TRUE, showWarnings = FALSE)

# ═════════════════════════════════════════
# STEP 1: Load & filter full transcriptome
# ═════════════════════════════════════════
message("=== STEP 1: Loading ===")
raw <- rio::import("data/raw/kidney_transcriptome.tsv")
n_genes_raw <- nrow(raw)
n_samples <- ncol(raw) - 1
message(sprintf("Raw: %d genes x %d samples", n_genes_raw, n_samples))

gene_symbols <- raw[[1]]
expr_data <- as.matrix(raw[, -1])
rownames(expr_data) <- gene_symbols
colnames(expr_data) <- colnames(raw)[-1]
storage.mode(expr_data) <- "numeric"

# Filtering report
n_ambiguous <- sum(grepl("\\?", gene_symbols))
expr_data <- expr_data[!grepl("\\?", gene_symbols), ]
n_after_ambiguous <- nrow(expr_data)

n_low <- sum(rowMeans(expr_data > 1) <= 0.1)
expr_data <- expr_data[rowMeans(expr_data > 1) > 0.1, ]
n_genes_final <- nrow(expr_data)

message(sprintf("Filtering: %d raw -> remove %d ambiguous symbols -> remove %d low expression = %d final",
                n_genes_raw, n_ambiguous, n_low, n_genes_final))

# ═════════════════════════════════════════
# STEP 2: Metadata & QC
# ═════════════════════════════════════════
message("\n=== STEP 2: Metadata ===")
samples <- colnames(expr_data)
meta <- tibble(
  sample = samples,
  study = ifelse(grepl("^GTEX-", samples), "GTEX", "TCGA"),
  sample_type = case_when(
    grepl("-01", samples, fixed = TRUE) ~ "Primary_Tumor",
    grepl("-11", samples, fixed = TRUE) ~ "Solid_Tissue_Normal",
    grepl("^GTEX-", samples) ~ "Normal_Tissue",
    TRUE ~ "Other"
  ),
  condition = case_when(
    sample_type == "Primary_Tumor" ~ "KIRP",
    study == "GTEX" ~ "Normal_GTEx",
    sample_type == "Solid_Tissue_Normal" ~ "TCGA_Normal",
    TRUE ~ "Other"
  ),
  participant = sub("^(TCGA-..-....).*", "\\1", samples)
)

kirp_participants <- unique(meta$participant[meta$condition == "KIRP"])
meta$paired_normal <- meta$condition == "TCGA_Normal" & meta$participant %in% kirp_participants
meta$normal_project <- ifelse(meta$paired_normal, "KIRP_adjacent", 
                               ifelse(meta$condition == "TCGA_Normal", "Other_RCC", NA))

message("Sample counts:")
print(table(meta$condition))
message(sprintf("TCGA Normals: %d (KIRP-adjacent: %d, Other projects: %d)", 
                sum(meta$condition == "TCGA_Normal"),
                sum(meta$paired_normal),
                sum(meta$condition == "TCGA_Normal" & !meta$paired_normal)))

# QC: PCA on full transcriptome (top 5000 most variable genes)
message("\n=== QC: Full transcriptome PCA ===")
vars <- apply(expr_data, 1, var)
top5k <- names(sort(vars, decreasing = TRUE))[1:min(5000, nrow(expr_data))]
pca <- prcomp(t(expr_data[top5k, ]), scale. = TRUE, center = TRUE)
pca_scores <- as.data.frame(pca$x)
pca_scores$condition <- meta$condition
pca_scores$study <- meta$study
var_pc1 <- round(summary(pca)$importance[2,1] * 100, 1)
var_pc2 <- round(summary(pca)$importance[2,2] * 100, 1)

png("results/v3/figures/PCA_transcriptome.png", width = 1000, height = 800, res = 130)
ggplot(pca_scores, aes(PC1, PC2, color = condition, shape = study)) +
  geom_point(alpha = 0.6, size = 2) +
  labs(title = "PCA — Full Transcriptome (top 5000 variable genes)",
       subtitle = sprintf("PC1: %.1f%% | PC2: %.1f%% | Colored by condition, shaped by study", var_pc1, var_pc2)) +
  theme_minimal(14)
dev.off()

# QC: Expression distributions
png("results/v3/figures/QC_density.png", width = 900, height = 600, res = 120)
plot(density(expr_data[, meta$condition == "KIRP"][, 1]), col = "#8A2BE2", lwd = 2,
     main = "Expression Density by Condition", xlab = "log2(expression)", ylim = c(0, 0.25))
lines(density(expr_data[, meta$condition == "TCGA_Normal"][, 1]), col = "#1F5BFF", lwd = 2)
lines(density(expr_data[, meta$condition == "Normal_GTEx"][, 1]), col = "#FFB347", lwd = 2)
legend("topright", c("KIRP", "TCGA Normal", "GTEx"), col = c("#8A2BE2", "#1F5BFF", "#FFB347"), lwd = 2)
dev.off()

# ═════════════════════════════════════════
# STEP 3: PRIMARY — Paired analysis
# ═════════════════════════════════════════
message("\n=== STEP 3: PRIMARY — Paired (32 pairs) ===")
paired_normals <- meta$sample[meta$paired_normal]
paired_tumors <- meta$sample[meta$condition == "KIRP" & meta$participant %in% meta$participant[meta$paired_normal]]
paired_meta <- meta[meta$sample %in% c(paired_normals, paired_tumors), ]
paired_meta$condition <- factor(paired_meta$condition, levels = c("TCGA_Normal", "KIRP"))
paired_meta$patient <- factor(paired_meta$participant)

E_paired <- expr_data[, paired_meta$sample, drop = FALSE]
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
message(sprintf("Paired DEGs: %d Up | %d Down", n_up, n_down))

# plotSA
png("results/v3/figures/plotSA_paired.png", width = 800, height = 600, res = 120)
plotSA(fit_paired, main = "Mean-Variance — Paired Analysis (31633 genes)")
dev.off()

# ═════════════════════════════════════════
# STEP 4: SECONDARY — KIRP vs KIRP-adjacent normals only (32 normals)
# ═════════════════════════════════════════
message("\n=== STEP 4: SECONDARY — 288 KIRP vs 32 adjacent normals ===")
tcga_kirp_idx <- meta$condition == "KIRP" | meta$paired_normal
tcga_kirp_meta <- meta[tcga_kirp_idx, ]
tcga_kirp_meta$condition <- factor(tcga_kirp_meta$condition, levels = c("TCGA_Normal", "KIRP"))
E_tcga_kirp <- expr_data[, tcga_kirp_meta$sample, drop = FALSE]

design_tcga_kirp <- model.matrix(~ condition, data = tcga_kirp_meta)
fit_tcga_kirp <- lmFit(E_tcga_kirp, design_tcga_kirp)
fit_tcga_kirp <- eBayes(fit_tcga_kirp, robust = TRUE, trend = TRUE)

deg_tcga_kirp <- topTable(fit_tcga_kirp, coef = "conditionKIRP", number = Inf, adjust.method = "BH", confint = TRUE)
deg_tcga_kirp$gene_id <- rownames(deg_tcga_kirp)
deg_tcga_kirp$regulation <- ifelse(deg_tcga_kirp$adj.P.Val < 0.05 & abs(deg_tcga_kirp$logFC) > 1,
                                    ifelse(deg_tcga_kirp$logFC > 0, "Up", "Down"), "NS")
message(sprintf("KIRP vs adjacent normals: %d Up | %d Down",
                sum(deg_tcga_kirp$regulation == "Up"), sum(deg_tcga_kirp$regulation == "Down")))

# ═════════════════════════════════════════
# STEP 5: EXPLORATORY — KIRP vs all TCGA normals (129)
# ═════════════════════════════════════════
message("\n=== STEP 5: EXPLORATORY — 288 KIRP vs 129 all TCGA normals ===")
tcga_all_idx <- meta$condition %in% c("KIRP", "TCGA_Normal")
tcga_all_meta <- meta[tcga_all_idx, ]
tcga_all_meta$condition <- factor(tcga_all_meta$condition, levels = c("TCGA_Normal", "KIRP"))
E_tcga_all <- expr_data[, tcga_all_meta$sample, drop = FALSE]

design_tcga_all <- model.matrix(~ condition, data = tcga_all_meta)
fit_tcga_all <- lmFit(E_tcga_all, design_tcga_all)
fit_tcga_all <- eBayes(fit_tcga_all, robust = TRUE, trend = TRUE)

deg_tcga_all <- topTable(fit_tcga_all, coef = "conditionKIRP", number = Inf, adjust.method = "BH", confint = TRUE)
deg_tcga_all$gene_id <- rownames(deg_tcga_all)
deg_tcga_all$regulation <- ifelse(deg_tcga_all$adj.P.Val < 0.05 & abs(deg_tcga_all$logFC) > 1,
                                   ifelse(deg_tcga_all$logFC > 0, "Up", "Down"), "NS")
message(sprintf("KIRP vs all TCGA normals: %d Up | %d Down",
                sum(deg_tcga_all$regulation == "Up"), sum(deg_tcga_all$regulation == "Down")))

# ═════════════════════════════════════════
# STEP 6: EXPLORATORY — KIRP vs GTEx (confounded)
# ═════════════════════════════════════════
message("\n=== STEP 6: EXPLORATORY — 288 KIRP vs 28 GTEx (confounded) ===")
gtx_idx <- meta$condition %in% c("KIRP", "Normal_GTEx")
gtx_meta <- meta[gtx_idx, ]
gtx_meta$condition <- factor(gtx_meta$condition, levels = c("Normal_GTEx", "KIRP"))
E_gtx <- expr_data[, gtx_meta$sample, drop = FALSE]

design_gtx <- model.matrix(~ condition, data = gtx_meta)
fit_gtx <- lmFit(E_gtx, design_gtx)
fit_gtx <- eBayes(fit_gtx, robust = TRUE, trend = TRUE)

deg_gtx <- topTable(fit_gtx, coef = "conditionKIRP", number = Inf, adjust.method = "BH", confint = TRUE)
deg_gtx$gene_id <- rownames(deg_gtx)
deg_gtx$regulation <- ifelse(deg_gtx$adj.P.Val < 0.05 & abs(deg_gtx$logFC) > 1,
                              ifelse(deg_gtx$logFC > 0, "Up", "Down"), "NS")
message(sprintf("KIRP vs GTEx: %d Up | %d Down (WARNING: confounded)",
                sum(deg_gtx$regulation == "Up"), sum(deg_gtx$regulation == "Down")))

# ═════════════════════════════════════════
# STEP 7: Paired plots for key genes
# ═════════════════════════════════════════
message("\n=== STEP 7: Paired plots ===")
key_genes <- c("ALDOB", "HK2", "PCK1", "G6PD", "TKT", "ADH1B", "FBP1")
for (g in intersect(key_genes, rownames(expr_data))) {
  normals <- paired_normals[paired_normals %in% colnames(expr_data)]
  tumors <- paired_tumors[paired_tumors %in% colnames(expr_data)]
  
  pd <- paired_meta
  pd$expression <- expr_data[g, pd$sample]
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
  hsa00010 = intersect(pw_membership$gene[pw_membership$in_hsa00010], rownames(expr_data)),
  hsa00030 = intersect(pw_membership$gene[pw_membership$in_hsa00030], rownames(expr_data)),
  hsa00020 = intersect(pw_membership$gene[pw_membership$in_hsa00020], rownames(expr_data))
)
names(pathways) <- c("hsa00010_Glycolysis", "hsa00030_PPP", "hsa00020_TCA")

# camera on paired analysis
camera_res <- camera(E_paired, pathways, design_paired, coef = coef_idx)
camera_res$pathway <- rownames(camera_res)
message("Camera gene set test (paired):")
print(camera_res)
rio::export(camera_res, "results/v3/tables/camera_gene_sets.csv")

# ═════════════════════════════════════════
# STEP 9: Extract pathway genes
# ═════════════════════════════════════════
message("\n=== STEP 9: Pathway extraction ===")
pw_names <- list(
  hsa00010 = "Glycolysis / Gluconeogenesis",
  hsa00030 = "Pentose Phosphate Pathway",
  hsa00020 = "Citrate Cycle (TCA)"
)

pw_results <- list()
for (pw_id in names(pw_names)) {
  pw_genes <- intersect(pw_membership$gene[pw_membership[[paste0("in_", pw_id)]]], rownames(expr_data))
  
  # Extract from all 4 analyses
  pw_paired <- deg_paired[deg_paired$gene_id %in% pw_genes, ]
  pw_tcga_kirp <- deg_tcga_kirp[deg_tcga_kirp$gene_id %in% pw_genes, ]
  pw_tcga_all <- deg_tcga_all[deg_tcga_all$gene_id %in% pw_genes, ]
  pw_gtx <- deg_gtx[deg_gtx$gene_id %in% pw_genes, ]
  
  pw_out <- data.frame(
    gene_id = sort(pw_genes),
    stringsAsFactors = FALSE
  )
  pw_out$logFC_Paired <- pw_paired$logFC[match(pw_out$gene_id, pw_paired$gene_id)]
  pw_out$CI.L_Paired <- pw_paired$CI.L[match(pw_out$gene_id, pw_paired$gene_id)]
  pw_out$CI.R_Paired <- pw_paired$CI.R[match(pw_out$gene_id, pw_paired$gene_id)]
  pw_out$FDR_Paired <- pw_paired$adj.P.Val[match(pw_out$gene_id, pw_paired$gene_id)]
  pw_out$AveExpr <- pw_paired$AveExpr[match(pw_out$gene_id, pw_paired$gene_id)]
  
  pw_out$logFC_TCGA_adj <- pw_tcga_kirp$logFC[match(pw_out$gene_id, pw_tcga_kirp$gene_id)]
  pw_out$FDR_TCGA_adj <- pw_tcga_kirp$adj.P.Val[match(pw_out$gene_id, pw_tcga_kirp$gene_id)]
  
  pw_out$logFC_TCGA_all <- pw_tcga_all$logFC[match(pw_out$gene_id, pw_tcga_all$gene_id)]
  pw_out$FDR_TCGA_all <- pw_tcga_all$adj.P.Val[match(pw_out$gene_id, pw_tcga_all$gene_id)]
  
  pw_out$logFC_GTEx <- pw_gtx$logFC[match(pw_out$gene_id, pw_gtx$gene_id)]
  pw_out$FDR_GTEx <- pw_gtx$adj.P.Val[match(pw_out$gene_id, pw_gtx$gene_id)]
  
  pw_out$pathway <- pw_id
  
  rio::export(pw_out, sprintf("results/v3/tables/DEG_%s.csv", pw_id))
  pw_results[[pw_id]] <- pw_out
  
  n_up <- sum(pw_out$FDR_Paired < 0.05 & pw_out$logFC_Paired > 1, na.rm = TRUE)
  n_down <- sum(pw_out$FDR_Paired < 0.05 & pw_out$logFC_Paired < -1, na.rm = TRUE)
  message(sprintf("  %s: %d genes, %d Up, %d Down", pw_id, nrow(pw_out), n_up, n_down))
}

# ═════════════════════════════════════════
# STEP 10: Supplementary table (all 106 genes)
# ═════════════════════════════════════════
message("\n=== STEP 10: Supplementary table ===")
all_pw_genes <- unique(unlist(lapply(pw_results, `[[`, "gene_id")))
supp_table <- pw_results[[1]][, c("gene_id", "pathway")]
for (pw_id in names(pw_names)[-1]) {
  other <- pw_results[[pw_id]][, c("gene_id", "pathway")]
  # Merge pathways for shared genes
  for (g in other$gene_id) {
    if (g %in% supp_table$gene_id) {
      supp_table$pathway[supp_table$gene_id == g] <- paste(supp_table$pathway[supp_table$gene_id == g], pw_id, sep = ";")
    } else {
      supp_table <- rbind(supp_table, other[other$gene_id == g, ])
    }
  }
}

# Merge all stats
supp <- supp_table
supp$logFC_Paired <- deg_paired$logFC[match(supp$gene_id, deg_paired$gene_id)]
supp$CI.L_Paired <- deg_paired$CI.L[match(supp$gene_id, deg_paired$gene_id)]
supp$CI.R_Paired <- deg_paired$CI.R[match(supp$gene_id, deg_paired$gene_id)]
supp$P.Value_Paired <- deg_paired$P.Value[match(supp$gene_id, deg_paired$gene_id)]
supp$FDR_Paired <- deg_paired$adj.P.Val[match(supp$gene_id, deg_paired$gene_id)]
supp$AveExpr <- deg_paired$AveExpr[match(supp$gene_id, deg_paired$gene_id)]
supp$logFC_TCGA_adj <- deg_tcga_kirp$logFC[match(supp$gene_id, deg_tcga_kirp$gene_id)]
supp$FDR_TCGA_adj <- deg_tcga_kirp$adj.P.Val[match(supp$gene_id, deg_tcga_kirp$gene_id)]
supp$logFC_TCGA_all <- deg_tcga_all$logFC[match(supp$gene_id, deg_tcga_all$gene_id)]
supp$FDR_TCGA_all <- deg_tcga_all$adj.P.Val[match(supp$gene_id, deg_tcga_all$gene_id)]
supp$logFC_GTEx <- deg_gtx$logFC[match(supp$gene_id, deg_gtx$gene_id)]
supp$FDR_GTEx <- deg_gtx$adj.P.Val[match(supp$gene_id, deg_gtx$gene_id)]

rio::export(supp, "results/v3/tables/Supplementary_Table_S1.csv")

# Unique DEGs
supp$is_DEG <- supp$FDR_Paired < 0.05 & abs(supp$logFC_Paired) > 1 & !is.na(supp$FDR_Paired)
n_unique_DEG <- sum(supp$is_DEG)
n_shared <- sum(supp$is_DEG & grepl(";", supp$pathway))
message(sprintf("Unique DEGs: %d (of which %d shared across pathways)", n_unique_DEG, n_shared))

# ═════════════════════════════════════════
# STEP 11: Concordance with proper metrics
# ═════════════════════════════════════════
message("\n=== STEP 11: Concordance ===")

calc_concordance <- function(x, y, label_x, label_y, gene_ids) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]; y <- y[ok]; ids <- gene_ids[ok]
  
  # Lin's CCC
  mx <- mean(x); my <- mean(y)
  vx <- var(x); vy <- var(y)
  sxy <- cov(x, y)
  ccc <- 2 * sxy / (vx + vy + (mx - my)^2)
  
  # Pearson
  r <- cor(x, y)
  
  # MAE, bias
  mae <- mean(abs(x - y))
  bias <- mean(x - y)
  
  # Regression
  lm_fit <- lm(y ~ x)
  slope <- coef(lm_fit)[2]
  intercept <- coef(lm_fit)[1]
  
  # Directional agreement (only genes with |logFC| > 0.5 in either)
  substantial <- abs(x) > 0.5 | abs(y) > 0.5
  dir_agree <- mean(sign(x[substantial]) == sign(y[substantial]))
  
  # Discordant genes table
  discordant_idx <- which(sign(x) != sign(y))
  n_discordant <- length(discordant_idx)
  
  list(ccc = ccc, r = r, mae = mae, bias = bias, slope = slope, intercept = intercept,
       dir_agree = dir_agree, n_total = length(x), n_discordant = n_discordant,
       label_x = label_x, label_y = label_y, discordant_genes = ids[discordant_idx])
}

# Paired vs TCGA-adjacent (secondary)
conc1 <- calc_concordance(supp$logFC_Paired, supp$logFC_TCGA_adj,
                           "Paired (32 pairs)", "KIRP vs adjacent normals (288 vs 32)",
                           supp$gene_id)
message(sprintf("\nPaired vs TCGA-adjacent:"))
message(sprintf("  Lin's CCC = %.3f | Pearson r = %.3f", conc1$ccc, conc1$r))
message(sprintf("  MAE = %.3f | Bias = %.3f | Slope = %.3f | Intercept = %.3f", 
                conc1$mae, conc1$bias, conc1$slope, conc1$intercept))
message(sprintf("  Directional (>|0.5|): %.1f%% | Discordant: %d", conc1$dir_agree * 100, conc1$n_discordant))

# Paired vs GTEx
conc2 <- calc_concordance(supp$logFC_Paired, supp$logFC_GTEx,
                           "Paired (32 pairs)", "KIRP vs GTEx (288 vs 28)",
                           supp$gene_id)
message(sprintf("\nPaired vs GTEx (confounded):"))
message(sprintf("  Lin's CCC = %.3f | Pearson r = %.3f", conc2$ccc, conc2$r))
message(sprintf("  MAE = %.3f | Bias = %.3f | Slope = %.3f | Intercept = %.3f",
                conc2$mae, conc2$bias, conc2$slope, conc2$intercept))
message(sprintf("  Directional (>|0.5|): %.1f%% | Discordant: %d", conc2$dir_agree * 100, conc2$n_discordant))

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
message(sprintf("Discordant genes saved: %d", nrow(discordant_df)))

# Bland-Altman plot
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

# Volcano plots
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
message(sprintf("Transcriptome: %d genes after filtering", n_genes_final))
message(sprintf("Primary: Paired analysis, 32 pairs, %d DEGs", sum(deg_paired$regulation != "NS")))
message(sprintf("Camera: %s", paste(rownames(camera_res), collapse=", ")))

# Save session info
writeLines(capture.output(sessionInfo()), "results/v3/sessionInfo.txt")
message("Session info saved")
