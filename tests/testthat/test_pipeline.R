library(testthat)

# Find repo root
repo_root <- normalizePath(file.path(dirname(getwd()), ".."), winslash = "/")
if (!dir.exists(file.path(repo_root, "data"))) {
  repo_root <- getwd()
}
message("Repo root: ", repo_root)

data_dir   <- file.path(repo_root, "data", "processed")
results_v3 <- file.path(repo_root, "results", "v3")
tables_v3  <- file.path(results_v3, "tables")

# ═══════════════════════════════════════
# PATHWAY MEMBERSHIP (shared v2/v3)
# ═══════════════════════════════════════
test_that("Pathway gene membership is correct", {
  pw <- rio::import(file.path(repo_root, "results", "tables", "pathway_gene_membership.csv"))
  expect_true(all(c("gene", "in_hsa00010", "in_hsa00030", "in_hsa00020", "in_matrix") %in% colnames(pw)))
  expect_equal(sum(pw$in_hsa00010), 67)
  expect_equal(sum(pw$in_hsa00030), 31)
  expect_equal(sum(pw$in_hsa00020), 30)
  expect_equal(sum(pw$in_matrix), 106)
})

# ═══════════════════════════════════════
# DEG FILES (v3)
# ═══════════════════════════════════════
test_that("DEG files exist and have correct structure", {
  for (pw_id in c("hsa00010", "hsa00020", "hsa00030")) {
    deg <- rio::import(file.path(tables_v3, paste0("DEG_", pw_id, ".csv")))
    expected_cols <- c("gene_id", "logFC_Paired", "FDR_Paired", "AveExpr",
                       "logFC_TCGA_adj", "logFC_GTEx", "pathway")
    expect_true(all(expected_cols %in% colnames(deg)))
  }
})

test_that("DEG counts match expected values", {
  deg10 <- rio::import(file.path(tables_v3, "DEG_hsa00010.csv"))
  deg20 <- rio::import(file.path(tables_v3, "DEG_hsa00020.csv"))
  deg30 <- rio::import(file.path(tables_v3, "DEG_hsa00030.csv"))
  
  expect_equal(nrow(deg10), 64)
  expect_equal(nrow(deg20), 29)
  expect_equal(nrow(deg30), 30)
  
  up10   <- sum(deg10$FDR_Paired < 0.05 & deg10$logFC_Paired > 1)
  down10 <- sum(deg10$FDR_Paired < 0.05 & deg10$logFC_Paired < -1)
  up20   <- sum(deg20$FDR_Paired < 0.05 & deg20$logFC_Paired > 1)
  down20 <- sum(deg20$FDR_Paired < 0.05 & deg20$logFC_Paired < -1)
  up30   <- sum(deg30$FDR_Paired < 0.05 & deg30$logFC_Paired > 1)
  down30 <- sum(deg30$FDR_Paired < 0.05 & deg30$logFC_Paired < -1)
  
  expect_equal(c(up10, down10), c(7, 19))
  expect_equal(c(up20, down20), c(0, 4))
  expect_equal(c(up30, down30), c(3, 9))
})

# ═══════════════════════════════════════
# KEY GENE VALUES
# ═══════════════════════════════════════
test_that("Key gene logFC values are consistent", {
  deg10 <- rio::import(file.path(tables_v3, "DEG_hsa00010.csv"))
  
  aldob <- deg10[deg10$gene_id == "ALDOB", ]
  hk2   <- deg10[deg10$gene_id == "HK2", ]
  pck1  <- deg10[deg10$gene_id == "PCK1", ]
  
  expect_true(aldob$logFC_Paired < -7)
  expect_true(hk2$logFC_Paired > 3)
  expect_true(pck1$logFC_Paired < -4)
  
  aldoa <- deg10[deg10$gene_id == "ALDOA", ]
  expect_true(aldoa$logFC_Paired > 1)
})

# ═══════════════════════════════════════
# CAMERA RESULTS
# ═══════════════════════════════════════
test_that("Camera results match expected values", {
  cam <- rio::import(file.path(tables_v3, "camera_gene_sets.csv"))
  expect_equal(nrow(cam), 3)
  
  tca_row <- cam[cam$pathway == "hsa00020_TCA", ]
  expect_true(tca_row$FDR < 0.01)
  expect_equal(tca_row$Direction, "Down")
  
  gly_row <- cam[cam$pathway == "hsa00010_Glycolysis", ]
  expect_true(gly_row$FDR > 0.05 & gly_row$FDR < 0.10)
})

# ═══════════════════════════════════════
# SUPPLEMENTARY TABLE S1
# ═══════════════════════════════════════
test_that("Supplementary Table S1 has correct structure", {
  s1 <- rio::import(file.path(tables_v3, "Supplementary_Table_S1.csv"))
  expect_equal(nrow(s1), 106)
  
  n_deg <- sum(s1$FDR_Paired < 0.05 & abs(s1$logFC_Paired) > 1, na.rm = TRUE)
  expect_equal(n_deg, 35)
  
  n_up   <- sum(s1$FDR_Paired < 0.05 & s1$logFC_Paired > 1, na.rm = TRUE)
  n_down <- sum(s1$FDR_Paired < 0.05 & s1$logFC_Paired < -1, na.rm = TRUE)
  expect_equal(c(n_up, n_down), c(9, 26))
})

# ═══════════════════════════════════════
# CONCORDANCE (discordant genes)
# ═══════════════════════════════════════
test_that("Discordant genes file has expected structure", {
  disc <- rio::import(file.path(tables_v3, "discordant_genes.csv"))
  expect_equal(nrow(disc), 45)
  expect_true(all(c("gene_id", "logFC_Paired", "logFC_GTEx") %in% colnames(disc)))
})

# ═══════════════════════════════════════
# FIGURES EXIST
# ═══════════════════════════════════════
test_that("Key figures exist", {
  fig_dir <- file.path(results_v3, "figures")
  expected_figs <- c("Volcano_hsa00010.png", "Volcano_hsa00020.png", 
                     "Volcano_hsa00030.png", "BlandAltman_Paired_vs_GTEx.png",
                     "PCA_transcriptome.png", "plotSA_paired.png")
  for (f in expected_figs) {
    expect_true(file.exists(file.path(fig_dir, f)))
  }
})

test_that("3D HTML files exist", {
  fig_dir <- file.path(results_v3, "figures")
  expect_true(file.exists(file.path(fig_dir, "PPI_network_3D.html")))
  expect_true(file.exists(file.path(fig_dir, "Volcano_hsa00010_3D.html")))
})

# ═══════════════════════════════════════
# REPOSITORY INTEGRITY
# ═══════════════════════════════════════
test_that("No v2 legacy files remain tracked", {
  skip_if_not(nchar(Sys.which("git")) > 0, "git not available in CI")
  v2_patterns <- c("results/differential_expression", "results/enrichment",
                   "results/figures/Volcano.png", "results/literature",
                   "scripts/01_data_provenance.R", "scripts/03_differential_expression.R")
  tracked <- system2("git", c("ls-files"), stdout = TRUE)
  for (pat in v2_patterns) {
    expect_false(any(grepl(pat, tracked, fixed = TRUE)))
  }
})

test_that("VERSION file is present and valid", {
  ver <- readLines(file.path(repo_root, "VERSION"))
  expect_true(length(ver) > 0)
  expect_match(ver[1], "v3\\.")
})
