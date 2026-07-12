library(testthat)

repo_root <- normalizePath(file.path(dirname(getwd()), ".."), winslash = "/")
if (!dir.exists(file.path(repo_root, "data"))) repo_root <- getwd()
message("Repo root: ", repo_root)

tables_v3  <- file.path(repo_root, "results", "v3", "tables")
figs_v3    <- file.path(repo_root, "results", "v3", "figures")

# Pathway membership
test_that("Pathway membership correct", {
  pw <- rio::import(file.path(repo_root, "results", "tables", "pathway_gene_membership.csv"))
  expect_equal(sum(pw$in_hsa00010), 67)
  expect_equal(sum(pw$in_hsa00030), 31)
  expect_equal(sum(pw$in_hsa00020), 30)
})

# DEG files
test_that("DEG files exist and counts match", {
  expect_equal(nrow(rio::import(file.path(tables_v3, "DEG_hsa00010.csv"))), 64)
  expect_equal(nrow(rio::import(file.path(tables_v3, "DEG_hsa00020.csv"))), 29)
  expect_equal(nrow(rio::import(file.path(tables_v3, "DEG_hsa00030.csv"))), 30)
})

# S1
test_that("Supplementary S1 correct", {
  s1 <- rio::import(file.path(tables_v3, "Supplementary_Table_S1.csv"))
  expect_equal(nrow(s1), 106)
  n_deg <- sum(s1$FDR_Paired < 0.05 & abs(s1$logFC_Paired) > 1, na.rm = TRUE)
  expect_equal(n_deg, 35)
})

# Camera
test_that("Camera results present", {
  cam <- rio::import(file.path(tables_v3, "camera_gene_sets.csv"))
  expect_equal(nrow(cam), 3)
  expect_true(all(c("hsa00020_TCA","hsa00030_PPP","hsa00010_Glycolysis") %in% cam$pathway))
})

# Enrichment tables exist
test_that("Enrichment tables exist", {
  for(f in c("ORA_KEGG_Up.csv","ORA_KEGG_Down.csv","ORA_Reactome_Up.csv","ORA_Reactome_Down.csv",
             "STRING_KEGG_Up.csv","STRING_KEGG_Down.csv","GSEA_KEGG.csv")) {
    expect_true(file.exists(file.path(tables_v3, f)))
  }
})

# Key figures exist
test_that("Key figures exist", {
  for(f in c("Volcano_hsa00010.png","Volcano_hsa00020.png","Volcano_hsa00030.png",
             "PPI_network.png","BlandAltman_Paired_vs_GTEx.png","Heatmap_35DEGs.png",
             "PCA_transcriptome.png","PPI_network_3D.html","Volcano_hsa00010_3D.html")) {
    expect_true(file.exists(file.path(figs_v3, f)))
  }
})

# Enrichment figures exist
test_that("Enrichment figures exist", {
  ef <- file.path(figs_v3, "enrichment")
  for(f in c("KEGG_Up_dotplot.png","KEGG_Down_dotplot.png","Reactome_Up_dotplot.png",
             "Reactome_Down_dotplot.png","STRING_Up_barplot.png","STRING_Down_barplot.png",
             "GSEA_KEGG_barplot.png")) {
    expect_true(file.exists(file.path(ef, f)))
  }
})

# VERSION
test_that("VERSION exists", {
  ver <- readLines(file.path(repo_root, "VERSION"))
  expect_match(ver[1], "v3\\.")
})

# kidney.tsv versioned
test_that("kidney.tsv is versioned", {
  expect_true(file.exists(file.path(repo_root, "data", "raw", "kidney.tsv")))
})
