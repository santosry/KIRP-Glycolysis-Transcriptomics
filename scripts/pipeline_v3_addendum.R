# pipeline_v3_addendum.R — análises faltantes
# Modelo não-pareado nas mesmas 64 amostras, treat(lfc=1), camera sensibilidade

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tibble)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
setwd(repo_root)

# ── Load data ──
message("Loading data...")
raw <- read.table("data/raw/kidney_transcriptome.tsv", header=TRUE, sep="\t", 
                   check.names=FALSE, row.names=1, comment.char="", quote="")
expr_data <- as.matrix(raw)
storage.mode(expr_data) <- "numeric"
expr_data <- expr_data[rowMeans(expr_data > 1) > 0.1, ]
message(sprintf("Genes: %d, Samples: %d", nrow(expr_data), ncol(expr_data)))

samples <- colnames(expr_data)
meta <- data.frame(
  sample = samples,
  condition = ifelse(grepl("-01", samples), "KIRP",
              ifelse(grepl("^GTEX-", samples), "Normal_GTEx", "TCGA_Normal")),
  participant = sub("^(TCGA-..-....).*", "\\1", samples),
  stringsAsFactors = FALSE
)

kirp_parts_all <- unique(meta$participant[meta$condition == "KIRP"])
norm_parts <- unique(meta$participant[meta$condition == "TCGA_Normal"])
paired_parts <- intersect(kirp_parts_all, norm_parts)
pn <- meta$sample[meta$condition == "TCGA_Normal" & meta$participant %in% paired_parts]
pt <- meta$sample[meta$condition == "KIRP" & meta$participant %in% paired_parts]

meta_64 <- meta[meta$sample %in% c(pn, pt), ]
meta_64$condition <- factor(meta_64$condition, levels = c("TCGA_Normal", "KIRP"))
meta_64$patient <- factor(meta_64$participant)
message(sprintf("Paired samples: %d", nrow(meta_64)))

E_64 <- expr_data[, meta_64$sample, drop = FALSE]

# ═══════════════════════════════
# 1. UNPAIRED MODEL ON SAME 64 SAMPLES
# ═══════════════════════════════
message("\n=== 1. UNPAIRED ON SAME 64 ===")
design_u64 <- model.matrix(~ condition, data = meta_64)
fit_u64 <- lmFit(E_64, design_u64)
fit_u64 <- eBayes(fit_u64, robust = TRUE, trend = TRUE)
deg_u64 <- topTable(fit_u64, coef = "conditionKIRP", number = Inf, adjust.method = "BH", confint = TRUE)
deg_u64$gene_id <- rownames(deg_u64)

# Compute paired model directly (no cached RDS)
dp <- model.matrix(~ patient + condition, data = meta_64)
fp <- lmFit(E_64, dp); fp <- eBayes(fp, robust = TRUE, trend = TRUE)
deg_paired <- topTable(fp, coef = grep("conditionKIRP", colnames(dp)), number = Inf, adjust.method = "BH", confint = TRUE)
deg_paired$gene_id <- rownames(deg_paired)

# Compare paired vs unpaired on same 64
common <- intersect(deg_paired$gene_id, deg_u64$gene_id)
x_p <- deg_paired$logFC[match(common, deg_paired$gene_id)]
y_u <- deg_u64$logFC[match(common, deg_u64$gene_id)]
ok <- !is.na(x_p) & !is.na(y_u); x_p <- x_p[ok]; y_u <- y_u[ok]

ccc_u64 <- 2*cov(x_p,y_u)/(var(x_p)+var(y_u)+(mean(x_p)-mean(y_u))^2)
mae_u64 <- mean(abs(x_p - y_u))
bias_u64 <- mean(x_p - y_u)  # paired - unpaired
r_u64 <- cor(x_p, y_u)
lmf <- lm(y_u ~ x_p)
mag_diff <- mean(abs(y_u) - abs(x_p))

message(sprintf("Paired vs Unpaired (same 64 samples, %d genes):", length(x_p)))
message(sprintf("  CCC = %.4f  r = %.4f  MAE = %.4f", ccc_u64, r_u64, mae_u64))
message(sprintf("  Bias (paired-unpaired) = %.4f", bias_u64))
message(sprintf("  Abs magnitude diff = %.4f", mag_diff))
message(sprintf("  Slope = %.4f  Intercept = %.4f", coef(lmf)[2], coef(lmf)[1]))

# Pathway genes only
pw <- read.csv("results/tables/pathway_gene_membership.csv", stringsAsFactors = FALSE)
pw_genes <- unique(pw$gene); pw_genes <- intersect(pw_genes, common)
xpw <- x_p[match(pw_genes, common)]; yuw <- y_u[match(pw_genes, common)]
okpw <- !is.na(xpw) & !is.na(yuw); xpw <- xpw[okpw]; yuw <- yuw[okpw]
ccc_pw_u64 <- 2*cov(xpw,yuw)/(var(xpw)+var(yuw)+(mean(xpw)-mean(yuw))^2)
subst <- abs(xpw) > 0.5 | abs(yuw) > 0.5
dir_pw <- mean(sign(xpw[subst]) == sign(yuw[subst]))
message(sprintf("  Pathway genes (%d): CCC=%.4f Dir(>0.5)=%d/%d=%.1f%%", 
                length(xpw), ccc_pw_u64, sum(sign(xpw[subst])==sign(yuw[subst])), 
                sum(subst), dir_pw*100))

# Save
unpaired64_results <- data.frame(
  metric = c("CCC","Pearson_r","MAE","Bias_signed","Abs_mag_diff","Slope","Intercept","N_genes","N_pathway","CCC_pathway","Dir_agree_pathway"),
  value = c(ccc_u64, r_u64, mae_u64, bias_u64, mag_diff, coef(lmf)[2], coef(lmf)[1], 
            length(x_p), length(xpw), ccc_pw_u64, dir_pw)
)
rio::export(unpaired64_results, "results/v3/tables/unpaired64_comparison.csv")

# ═══════════════════════════════
# 2. treat(lfc=1) SENSITIVITY
# ═══════════════════════════════
message("\n=== 2. treat(lfc=1) ===")
# treat works on the eBayes fit; use the last coefficient (conditionKIRP)
coef_idx <- ncol(dp)
treat_res <- treat(fp, lfc = 1)
treat_deg <- topTreat(treat_res, number = Inf, adjust.method = "BH")
treat_deg$gene_id <- rownames(treat_deg)

# Compare with standard DEG classification
standard_deg <- deg_paired
standard_deg$is_DEG <- standard_deg$adj.P.Val < 0.05 & abs(standard_deg$logFC) > 1
treat_deg$is_DEG_treat <- treat_deg$adj.P.Val < 0.05

common2 <- intersect(standard_deg$gene_id, treat_deg$gene_id)
comp <- data.frame(
  gene_id = common2,
  standard_DEG = standard_deg$is_DEG[match(common2, standard_deg$gene_id)],
  treat_DEG = treat_deg$is_DEG_treat[match(common2, treat_deg$gene_id)],
  stringsAsFactors = FALSE
)
n_std <- sum(comp$standard_DEG)
n_treat <- sum(comp$treat_DEG)
n_both <- sum(comp$standard_DEG & comp$treat_DEG)
n_std_only <- sum(comp$standard_DEG & !comp$treat_DEG)
n_treat_only <- sum(!comp$standard_DEG & comp$treat_DEG)

message(sprintf("Standard DEGs (FDR<0.05 & |logFC|>1): %d", n_std))
message(sprintf("treat DEGs (FDR<0.05 for |logFC|>1): %d", n_treat))
message(sprintf("Both: %d | Standard only: %d | treat only: %d", n_both, n_std_only, n_treat_only))

# Pathway genes
pw_treat <- comp[comp$gene_id %in% pw_genes, ]
n_pw_std <- sum(pw_treat$standard_DEG)
n_pw_treat <- sum(pw_treat$treat_DEG)
n_pw_both <- sum(pw_treat$standard_DEG & pw_treat$treat_DEG)
message(sprintf("Pathway genes: Std DEG=%d, treat DEG=%d, both=%d", n_pw_std, n_pw_treat, n_pw_both))

rio::export(comp, "results/v3/tables/treat_sensitivity.csv")

# ═══════════════════════════════
# 3. camera WITH inter.gene.cor=NA
# ═══════════════════════════════
message("\n=== 3. camera with inter.gene.cor=NA ===")
pathways <- list(
  hsa00010_Glycolysis = intersect(pw$gene[pw$in_hsa00010], rownames(E_64)),
  hsa00030_PPP = intersect(pw$gene[pw$in_hsa00030], rownames(E_64)),
  hsa00020_TCA = intersect(pw$gene[pw$in_hsa00020], rownames(E_64))
)

cam_na <- camera(E_64, pathways, dp, coef = grep("conditionKIRP", colnames(dp)), 
                 inter.gene.cor = NA)
cam_na <- as.data.frame(cam_na)
cam_na$pathway <- rownames(cam_na)

# Compare with inter.gene.cor=0.01
cam_01 <- read.csv("results/v3/tables/camera_gene_sets.csv", stringsAsFactors = FALSE)

message("Camera comparison (0.01 vs NA):")
for (i in 1:nrow(cam_na)) {
  pw_name <- cam_na$pathway[i]
  fdr_01 <- cam_01$FDR[cam_01$pathway == pw_name]
  fdr_na <- cam_na$FDR[i]
  message(sprintf("  %s: FDR(0.01)=%.4f FDR(NA)=%.4f", pw_name, fdr_01, fdr_na))
}

rio::export(cam_na, "results/v3/tables/camera_inter_gene_cor_NA.csv")

# ═══════════════════════════════
# 4. SESSION INFO & MANIFEST
# ═══════════════════════════════
message("\n=== 4. Session info ===")
writeLines(capture.output(sessionInfo()), "results/v3/sessionInfo_full.txt")

# Generate manifest
files <- list.files("results/v3", recursive = TRUE, full.names = TRUE)
manifest <- data.frame(
  file = files,
  size_bytes = file.info(files)$size,
  last_modified = file.info(files)$mtime,
  stringsAsFactors = FALSE
)
rio::export(manifest, "results/v3/results_manifest.csv")

# Generate checksums
if (requireNamespace("digest", quietly = TRUE)) {
  hashes <- sapply(files, function(f) digest::digest(f, algo = "sha256", file = TRUE))
  writeLines(paste(hashes, basename(files)), "results/v3/checksums_sha256.txt")
}

message("\n=== ADDENDUM COMPLETE ===")
message(sprintf("Unpaired-64 CCC: %.4f", ccc_u64))
message(sprintf("treat pathway DEGs: %d (vs standard %d)", n_pw_treat, n_pw_std))
