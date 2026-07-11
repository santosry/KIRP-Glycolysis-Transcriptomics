# 15_tcga_only_sensitivity.R
# TCGA-ONLY SENSITIVITY ANALYSIS: KIRP vs TCGA_Normal (adjacent normal)
# This avoids the cohort confounding present in the main KIRP vs GTEx comparison

suppressPackageStartupMessages({
  library(limma); library(dplyr); library(tibble); library(rio); library(ggplot2)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)

tcga_meta_file <- file.path(repo_root, "data", "processed", "metadata_tcga_only.rds")
tcga_expr_file <- file.path(repo_root, "data", "processed", "expression_matrix_tcga_only.rds")

if (!file.exists(tcga_meta_file) || !file.exists(tcga_expr_file)) {
  stop("TCGA-only data not found. Run 02_prepare_data.R first.")
}

dir.create(file.path(repo_root, "results", "sensitivity"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)

meta_tcga <- readRDS(tcga_meta_file)
expr_tcga <- readRDS(tcga_expr_file)

# Check for missing values
if (anyNA(expr_tcga)) {
  message("Missing values in TCGA expression matrix - will use complete.cases")
}

message("TCGA-only analysis: ", nrow(meta_tcga), " samples")
message("  KIRP: ", sum(meta_tcga$proposed_condition == "KIRP"))
message("  TCGA Normal: ", sum(meta_tcga$proposed_condition == "TCGA_Normal"))

# ── Check TCGA Normal provenance ──
# Inspect which TCGA projects the normals come from
tcga_normals <- meta_tcga |> filter(proposed_condition == "TCGA_Normal")
message("\nTCGA Normal samples by study:")
print(table(tcga_normals$study))

# ── Expression matrix ──
E_tcga <- t(expr_tcga)
storage.mode(E_tcga) <- "numeric"

# Remove zero-variance genes
gene_vars <- apply(E_tcga, 1, var, na.rm = TRUE)
E_tcga <- E_tcga[gene_vars > 0, , drop = FALSE]
n_genes <- nrow(E_tcga)
message("\nGenes with variance > 0: ", n_genes)

# ── Design ──
meta_tcga$condition <- factor(meta_tcga$proposed_condition, 
                               levels = c("TCGA_Normal", "KIRP"))

# Check study confounding within TCGA
if (length(unique(meta_tcga$study)) > 1) {
  message("\nTCGA-only study distribution:")
  print(table(meta_tcga$condition, meta_tcga$study))
  # For TCGA-only, study should be all "TCGA"
  design_tcga <- model.matrix(~ 0 + condition, data = meta_tcga)
} else {
  design_tcga <- model.matrix(~ 0 + condition, data = meta_tcga)
}

colnames(design_tcga) <- gsub("condition", "", colnames(design_tcga))

# ── limma ──
message("\nFitting limma model (TCGA-only)...")
fit_tcga <- lmFit(E_tcga, design_tcga)
contrast_tcga <- makeContrasts(KIRP_vs_TCGA_Normal = KIRP - TCGA_Normal, levels = design_tcga)
fit2_tcga <- contrasts.fit(fit_tcga, contrast_tcga)
fit2_tcga <- eBayes(fit2_tcga, robust = TRUE, trend = TRUE)

deg_tcga <- topTable(fit2_tcga, number = Inf, adjust.method = "BH", sort.by = "P")
deg_tcga <- as.data.frame(deg_tcga) |> 
  tibble::rownames_to_column("gene_id") |> 
  as_tibble()

# Classify
deg_tcga$regulation <- "NS"
deg_tcga$regulation[deg_tcga$adj.P.Val < 0.05 & deg_tcga$logFC > 1] <- "Up"
deg_tcga$regulation[deg_tcga$adj.P.Val < 0.05 & deg_tcga$logFC < -1] <- "Down"

n_up_tcga <- sum(deg_tcga$regulation == "Up")
n_down_tcga <- sum(deg_tcga$regulation == "Down")
message("TCGA-only DEGs (|logFC|>1, FDR<0.05): ", n_up_tcga, " Up | ", n_down_tcga, " Down")

# ── Compare with main analysis ──
deg_main_file <- file.path(repo_root, "results", "differential_expression", "DEG_global.csv")
if (file.exists(deg_main_file)) {
  deg_main <- rio::import(deg_main_file)
  
  comparison <- deg_main |> 
    select(gene_id, logFC_main = logFC, adj.P.Val_main = adj.P.Val, regulation_main = regulation) |>
    full_join(deg_tcga |> select(gene_id, logFC_tcga = logFC, adj.P.Val_tcga = adj.P.Val, regulation_tcga = regulation),
              by = "gene_id")
  
  # logFC correlation
  common <- comparison |> filter(!is.na(logFC_main) & !is.na(logFC_tcga))
  if (nrow(common) > 0) {
    lfc_cor <- cor(common$logFC_main, common$logFC_tcga)
    dir_agree <- mean(sign(common$logFC_main) == sign(common$logFC_tcga))
    message(sprintf("\nMain vs TCGA-only: logFC correlation = %.4f | Directional agreement = %.1f%%",
                    lfc_cor, 100 * dir_agree))
  }
  
  rio::export(comparison, file.path(repo_root, "results", "sensitivity", "main_vs_tcga_only.csv"))
  
  # ── Volcano comparison ──
  deg_tcga$dataset <- "TCGA-only"
  deg_main$dataset <- "Main (KIRP vs GTEx)"
  
  deg_tcga$color <- deg_tcga$regulation
  deg_main$color <- deg_main$regulation
  
  p_comp <- bind_rows(
    deg_main |> select(gene_id, logFC, adj.P.Val, regulation, dataset),
    deg_tcga |> select(gene_id, logFC, adj.P.Val, regulation, dataset)
  ) |>
    mutate(color = regulation) |>
    ggplot(aes(logFC, -log10(adj.P.Val), color = color)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_color_manual(values = c(Up = "#1F5BFF", Down = "#8A2BE2", NS = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3) +
    facet_wrap(~ dataset, ncol = 1) +
    labs(title = "Sensitivity: Main vs TCGA-only",
         subtitle = paste0("Main: ", sum(deg_main$regulation != "NS"), " DEGs | TCGA-only: ", n_up_tcga + n_down_tcga, " DEGs"),
         x = "log2 Fold Change", y = "-log10(FDR)") +
    theme_classic(14)
  ggsave(file.path(repo_root, "results", "figures", "Volcano_main_vs_tcga.png"),
         p_comp, width = 9, height = 10, dpi = 300)
  
  # Save combined comparison table
  top_genes <- c("ALDOB", "HK2", "ALDOA", "PKM", "PFKP", "GAPDH", "ENO1", "PCK1",
                 "ADH1B", "ADH1C", "ALDH3B1", "ALDH3B2", "ADH1A", "ADH4", "ADH6")
  comparison_table <- comparison |>
    filter(gene_id %in% top_genes) |>
    arrange(desc(abs(logFC_main)))
  rio::export(comparison_table, file.path(repo_root, "results", "sensitivity", "tcga_vs_main_top_genes.csv"))
  message("\nTop genes comparison:")
  for (i in 1:nrow(comparison_table)) {
    r <- comparison_table[i, ]
    message(sprintf("  %-10s Main: % 6.2f (%s) | TCGA: % 6.2f (%s)",
                    r$gene_id, r$logFC_main, r$regulation_main, r$logFC_tcga, r$regulation_tcga))
  }
}

# Save TCGA-only DEG results
rio::export(deg_tcga, file.path(repo_root, "results", "sensitivity", "DEG_TCGA_only.csv"))

message("\n✓ TCGA-only sensitivity analysis complete.")
