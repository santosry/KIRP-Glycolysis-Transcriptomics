library(testthat)

# Find repo root
repo_root <- normalizePath(file.path(dirname(getwd()), ".."), winslash = "/")
if (!dir.exists(file.path(repo_root, "data"))) {
  # Try as child directory
  repo_root <- getwd()
}
message("Repo root: ", repo_root)

data_dir <- file.path(repo_root, "data", "processed")
results_dir <- file.path(repo_root, "results")

test_that("Expression matrix dimensions are correct", {
  expr <- readRDS(file.path(data_dir, "expression_matrix.rds"))
  meta <- readRDS(file.path(data_dir, "metadata.rds"))
  
  expect_equal(nrow(expr), 316)
  expect_equal(nrow(meta), 316)
  expect_equal(ncol(expr), 66)
  expect_equal(nrow(expr), nrow(meta))
})

test_that("No missing values in expression matrix", {
  expr <- readRDS(file.path(data_dir, "expression_matrix.rds"))
  expect_equal(sum(is.na(expr)), 0)
})

test_that("Sample classification is correct", {
  meta <- readRDS(file.path(data_dir, "metadata.rds"))
  
  expect_equal(sum(meta$condition == "KIRP"), 288)
  expect_equal(sum(meta$condition == "Normal_GTEx"), 28)
  expect_setequal(levels(meta$condition), c("Normal_GTEx", "KIRP"))
})

test_that("DEG results have correct structure", {
  deg <- rio::import(file.path(results_dir, "differential_expression", "DEG_global.csv"))
  
  expect_true(all(c("gene_id", "logFC", "adj.P.Val", "regulation", "AveExpr", "t") %in% colnames(deg)))
  expect_equal(nrow(deg), 66)
  expect_setequal(unique(deg$regulation), c("NS", "Up", "Down"))
  
  expect_equal(sum(deg$regulation == "Up"), 19)
  expect_equal(sum(deg$regulation == "Down"), 12)
  expect_equal(sum(deg$regulation == "NS"), 35)
})

test_that("Key gene logFC values are consistent", {
  deg <- rio::import(file.path(results_dir, "differential_expression", "DEG_global.csv"))
  
  aldob <- deg[deg$gene_id == "ALDOB", ]
  hk2 <- deg[deg$gene_id == "HK2", ]
  
  expect_true(aldob$logFC < -7)
  expect_true(hk2$logFC > 3)
  expect_equal(aldob$regulation, "Down")
  expect_equal(hk2$regulation, "Up")
  
  # Check ALDOA and ALDH divergences
  aldoa <- deg[deg$gene_id == "ALDOA", ]
  expect_gt(aldoa$logFC, 1)
  
  aldh3b1 <- deg[deg$gene_id == "ALDH3B1", ]
  aldh3b2 <- deg[deg$gene_id == "ALDH3B2", ]
  expect_gt(aldh3b1$logFC, 0)
  expect_lt(aldh3b2$logFC, 0)
})

test_that("Design matrix has correct confounding structure", {
  meta <- readRDS(file.path(data_dir, "metadata.rds"))
  
  kirp_studies <- unique(meta$study[meta$condition == "KIRP"])
  normal_studies <- unique(meta$study[meta$condition == "Normal_GTEx"])
  
  expect_length(intersect(kirp_studies, normal_studies), 0)
  
  design <- model.matrix(~ condition + study, data = meta)
  expect_equal(qr(design)$rank, 2)
  expect_equal(ncol(design), 3)
})

test_that("Scale is log2 (not linear)", {
  expr <- readRDS(file.path(data_dir, "expression_matrix.rds"))
  vals <- as.vector(expr)
  vals <- vals[is.finite(vals)]
  
  expect_lt(max(vals), 25)
  expect_gt(median(vals), 9)
})

test_that("No duplicate sample IDs", {
  meta <- readRDS(file.path(data_dir, "metadata.rds"))
  expect_false(any(duplicated(meta$sample)))
})

test_that("TCGA-only dataset exists and has correct structure", {
  skip_if_not(file.exists(file.path(data_dir, "metadata_tcga_only.rds")))
  meta_t <- readRDS(file.path(data_dir, "metadata_tcga_only.rds"))
  expr_t <- readRDS(file.path(data_dir, "expression_matrix_tcga_only.rds"))
  
  expect_equal(nrow(meta_t), 417)
  expect_equal(sum(meta_t$proposed_condition == "KIRP"), 288)
  expect_equal(sum(meta_t$proposed_condition == "TCGA_Normal"), 129)
})
