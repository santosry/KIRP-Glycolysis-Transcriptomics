# 04_v2_main_analysis.R
# v2.0.0 CENTRAL CARBON METABOLISM ANALYSIS
# Comparator hierarchy: TCGA-paired > TCGA-all > GTEx exploratory

suppressPackageStartupMessages({
  library(limma); library(dplyr); library(tibble); library(rio)
  library(ggplot2); library(ggrepel); library(pheatmap); library(RColorBrewer)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)

# === LOAD DATA ===
raw_file <- file.path(repo_root, "data", "raw", "kidney.tsv")
kidney <- rio::import(raw_file)

# === COLUMN HANDLING (definitive) ===
meta_cols <- c("sample","samples","TCGA_GTEX_main_category","_sample_type","_study","_primary_site","OS","OS.time")
# Rename underscore-prefixed
colnames(kidney)[colnames(kidney) == "_sample_type"] <- "sample_type"
colnames(kidney)[colnames(kidney) == "_study"] <- "study"
colnames(kidney)[colnames(kidney) == "_primary_site"] <- "primary_site"
gene_cols <- setdiff(colnames(kidney), c("sample","samples","TCGA_GTEX_main_category","sample_type","study","primary_site","OS","OS.time"))

# === EXPRESSION MATRIX ===
expr <- as.matrix(kidney[, gene_cols, drop = FALSE])
storage.mode(expr) <- "numeric"
# Already log2(norm_count+1) — verified
message("Expression matrix: ", nrow(expr), " samples x ", ncol(expr), " genes")
message("Scale: log2(norm_count+1) | Range: [", min(expr), ", ", max(expr), "]")

# === COMPARATOR HIERARCHY ===
# 1. KIRP tumor (primary tumor, TCGA KIRP category)
# 2. GTEx normal kidney
# 3. TCGA solid tissue normal (all 129, from various renal projects)
# 4. TCGA normal matched to KIRP participants (32, potential paired)

meta <- kidney[, c("sample", "TCGA_GTEX_main_category", "sample_type", "study")]

meta$group <- NA
meta$group[meta$TCGA_GTEX_main_category == "TCGA Kidney Papillary Cell Carcinoma"] <- "KIRP_Tumor"
meta$group[meta$TCGA_GTEX_main_category == "GTEX Kidney"] <- "GTEx_Normal"
meta$group[meta$study == "TCGA" & meta$sample_type == "Solid Tissue Normal" & 
            (is.na(meta$TCGA_GTEX_main_category) | meta$TCGA_GTEX_main_category == "")] <- "TCGA_Normal"

# Identify KIRP-matched normals
kirp_participants <- meta |>
  filter(group == "KIRP_Tumor") |>
  mutate(participant = sapply(strsplit(sample, "-"), function(x) x[3])) |>
  pull(participant)

meta$participant <- sapply(strsplit(meta$sample, "-"), function(x) x[3])
meta$matched_kirp <- meta$participant %in% kirp_participants
meta$group[meta$group == "TCGA_Normal" & meta$matched_kirp] <- "TCGA_Normal_KIRP"
meta$group[meta$group == "TCGA_Normal" & !meta$matched_kirp] <- "TCGA_Normal_Other"

# Filter to classified
meta <- meta[!is.na(meta$group), ]
expr <- expr[match(meta$sample, kidney$sample), , drop = FALSE]

# === COUNTS ===
counts <- table(meta$group)
message("\n=== SAMPLE GROUPS ===")
for (nm in names(counts)) message(sprintf("  %-25s: %3d", nm, counts[nm]))

# === PAIRED ANALYSIS (gold standard) ===
# 32 KIRP participants have both tumor and KIRP-matched normal
paired_participants <- meta |>
  filter(group %in% c("KIRP_Tumor", "TCGA_Normal_KIRP")) |>
  group_by(participant) |>
  filter(n() == 2) |>
  pull(participant) |> unique()

message("\n=== PAIRED ANALYSIS ===")
message("Participants with both tumor and normal: ", length(paired_participants))

paired_meta <- meta |> filter(participant %in% paired_participants)
paired_expr <- expr[meta$participant %in% paired_participants, , drop = FALSE]

# Build paired design
paired_meta$patient <- factor(paired_meta$participant)
paired_meta$condition <- factor(ifelse(paired_meta$group == "KIRP_Tumor", "Tumor", "Normal"),
                                 levels = c("Normal", "Tumor"))

design_paired <- model.matrix(~ patient + condition, data = paired_meta)

E_paired <- t(paired_expr)
fit_paired <- lmFit(E_paired, design_paired)
fit_paired <- eBayes(fit_paired, robust = TRUE, trend = TRUE)
deg_paired <- topTable(fit_paired, coef = "conditionTumor", number = Inf, adjust.method = "BH", sort.by = "P")
deg_paired <- as.data.frame(deg_paired) |> tibble::rownames_to_column("gene_id") |> as_tibble()
deg_paired$regulation <- "NS"
deg_paired$regulation[deg_paired$adj.P.Val < 0.05 & deg_paired$logFC > 1] <- "Up"
deg_paired$regulation[deg_paired$adj.P.Val < 0.05 & deg_paired$logFC < -1] <- "Down"

message("Paired DEGs: ", sum(deg_paired$regulation == "Up"), " Up | ", 
        sum(deg_paired$regulation == "Down"), " Down | ",
        sum(deg_paired$regulation != "NS"), " total")

# === UNPAIRED ANALYSIS — Comparator 1: KIRP vs GTEx (cross-cohort, exploratory) ===
meta_gtx <- meta |> filter(group %in% c("KIRP_Tumor", "GTEx_Normal"))
meta_gtx$cond <- factor(meta_gtx$group, levels = c("GTEx_Normal", "KIRP_Tumor"))
expr_gtx <- expr[match(meta_gtx$sample, meta$sample), , drop = FALSE]

E_gtx <- t(expr_gtx)
design_gtx <- model.matrix(~ cond, data = meta_gtx)
fit_gtx <- lmFit(E_gtx, design_gtx)
fit_gtx <- eBayes(fit_gtx, robust = TRUE, trend = TRUE)
deg_gtx <- topTable(fit_gtx, coef = 2, number = Inf, adjust.method = "BH", sort.by = "P")
deg_gtx <- as.data.frame(deg_gtx) |> tibble::rownames_to_column("gene_id") |> as_tibble()
deg_gtx$regulation <- "NS"
deg_gtx$regulation[deg_gtx$adj.P.Val < 0.05 & deg_gtx$logFC > 1] <- "Up"
deg_gtx$regulation[deg_gtx$adj.P.Val < 0.05 & deg_gtx$logFC < -1] <- "Down"
message("KIRP vs GTEx DEGs: ", sum(deg_gtx$regulation == "Up"), " Up | ",
        sum(deg_gtx$regulation == "Down"), " Down")

# === UNPAIRED ANALYSIS — Comparator 2: KIRP vs all TCGA Normal ===
meta_tcga <- meta |> filter(group %in% c("KIRP_Tumor", "TCGA_Normal_KIRP", "TCGA_Normal_Other"))
meta_tcga$cond <- factor(ifelse(meta_tcga$group == "KIRP_Tumor", "Tumor", "Normal"),
                          levels = c("Normal", "Tumor"))
expr_tcga <- expr[match(meta_tcga$sample, meta$sample), , drop = FALSE]

E_tcga <- t(expr_tcga)
design_tcga <- model.matrix(~ cond, data = meta_tcga)
fit_tcga <- lmFit(E_tcga, design_tcga)
fit_tcga <- eBayes(fit_tcga, robust = TRUE, trend = TRUE)
deg_tcga <- topTable(fit_tcga, coef = 2, number = Inf, adjust.method = "BH", sort.by = "P")
deg_tcga <- as.data.frame(deg_tcga) |> tibble::rownames_to_column("gene_id") |> as_tibble()
deg_tcga$regulation <- "NS"
deg_tcga$regulation[deg_tcga$adj.P.Val < 0.05 & deg_tcga$logFC > 1] <- "Up"
deg_tcga$regulation[deg_tcga$adj.P.Val < 0.05 & deg_tcga$logFC < -1] <- "Down"
message("KIRP vs all TCGA Normal DEGs: ", sum(deg_tcga$regulation == "Up"), " Up | ",
        sum(deg_tcga$regulation == "Down"), " Down")

# === ROBUSTNESS CLASSIFICATION ===
# Pre-registered criteria:
# ROBUST: FDR<0.05 & |logFC|>1 in ALL three analyses (paired, GTEx, TCGA) AND same direction
# PARCIALLY ROBUST: FDR<0.05 & |logFC|>1 in at least 2 of 3, same direction
# SENSITIVE: FDR<0.05 & |logFC|>1 in only 1 of 3, or direction disagreement
# INCONCLUSIVE: FDR>=0.05 in all three

robustness <- tibble(
  gene_id = gene_cols
) |>
  left_join(deg_paired |> select(gene_id, logFC_paired = logFC, FDR_paired = adj.P.Val, reg_paired = regulation), by = "gene_id") |>
  left_join(deg_gtx |> select(gene_id, logFC_gtx = logFC, FDR_gtx = adj.P.Val, reg_gtx = regulation), by = "gene_id") |>
  left_join(deg_tcga |> select(gene_id, logFC_tcga = logFC, FDR_tcga = adj.P.Val, reg_tcga = regulation), by = "gene_id")

robustness <- robustness |>
  mutate(
    sig_paired = !is.na(FDR_paired) & FDR_paired < 0.05 & abs(logFC_paired) > 1,
    sig_gtx    = !is.na(FDR_gtx) & FDR_gtx < 0.05 & abs(logFC_gtx) > 1,
    sig_tcga   = !is.na(FDR_tcga) & FDR_tcga < 0.05 & abs(logFC_tcga) > 1,
    n_sig = sig_paired + sig_gtx + sig_tcga,
    dir_paired = sign(logFC_paired),
    dir_gtx = sign(logFC_gtx),
    dir_tcga = sign(logFC_tcga),
    dir_agree = (dir_paired == dir_gtx | is.na(dir_paired) | is.na(dir_gtx)) &
                (dir_paired == dir_tcga | is.na(dir_paired) | is.na(dir_tcga)) &
                (dir_gtx == dir_tcga | is.na(dir_gtx) | is.na(dir_tcga)),
    robustness_class = case_when(
      n_sig == 3 & dir_agree ~ "ROBUST",
      n_sig >= 2 & dir_agree ~ "PARTIALLY_ROBUST",
      n_sig >= 1 & !dir_agree ~ "DIRECTION_DISCORDANT",
      n_sig >= 1 & dir_agree ~ "SENSITIVE_TO_COMPARATOR",
      TRUE ~ "INCONCLUSIVE"
    )
  )

# Count by class
rob_counts <- table(robustness$robustness_class)
message("\n=== ROBUSTNESS CLASSIFICATION ===")
for (nm in names(rob_counts)) message(sprintf("  %-25s: %2d", nm, rob_counts[nm]))

# Save all results
dir.create(file.path(repo_root, "results", "v2"), recursive = TRUE, showWarnings = FALSE)
rio::export(deg_paired, file.path(repo_root, "results", "v2", "DEG_paired.csv"))
rio::export(deg_gtx, file.path(repo_root, "results", "v2", "DEG_gtx.csv"))
rio::export(deg_tcga, file.path(repo_root, "results", "v2", "DEG_tcga_all.csv"))
rio::export(robustness, file.path(repo_root, "results", "v2", "robustness_classification.csv"))
rio::export(meta, file.path(repo_root, "results", "v2", "sample_groups.csv"))

# === PRINT KEY GENES ===
message("\n=== ROBUST GENES ===")
rb <- robustness |> filter(robustness_class == "ROBUST") |> arrange(desc(abs(logFC_gtx)))
for (i in 1:nrow(rb)) {
  message(sprintf("  %-10s Paired:% 6.2f GTEx:% 6.2f TCGA:% 6.2f",
                  rb$gene_id[i], rb$logFC_paired[i], rb$logFC_gtx[i], rb$logFC_tcga[i]))
}

message("\n✓ v2 main analysis complete.")
