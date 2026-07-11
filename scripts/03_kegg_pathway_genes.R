# Get all three KEGG pathways
library(KEGGREST)
library(dplyr)
library(rio)

get_kegg_genes <- function(pathway) {
  pw <- keggGet(pathway)
  genes_raw <- pw[[1]]$GENE
  symbols <- genes_raw[seq(2, length(genes_raw), by = 2)]
  symbols <- gsub("\\s*\\[.*$", "", symbols)
  symbols <- trimws(symbols)
  symbols <- gsub(";.*$", "", symbols)
  symbols <- trimws(symbols)
  symbols <- unique(symbols)
  symbols <- sort(symbols)
  symbols
}

hsa00010 <- get_kegg_genes("hsa00010")
hsa00030 <- get_kegg_genes("hsa00030")
hsa00020 <- get_kegg_genes("hsa00020")

cat("hsa00010 (Glycolysis/Gluconeogenesis):", length(hsa00010), "genes\n")
cat("hsa00030 (Pentose Phosphate Pathway):  ", length(hsa00030), "genes\n")
cat("hsa00020 (Citrate Cycle / TCA):         ", length(hsa00020), "genes\n")

all_three <- unique(c(hsa00010, hsa00030, hsa00020))
cat("Union of all three pathways:", length(all_three), "unique genes\n")

# What's in kidney.tsv?
raw_file <- "data/raw/kidney.tsv"
hdr <- strsplit(readLines(raw_file, 1), "\t")[[1]]
meta_cols <- c("sample","samples","TCGA_GTEX_main_category","_sample_type","_study","_primary_site","OS","OS.time")
matrix_genes <- setdiff(hdr, meta_cols)
cat("Genes in kidney.tsv:", length(matrix_genes), "\n")

in_00010 <- intersect(hsa00010, matrix_genes)
in_00030 <- intersect(hsa00030, matrix_genes)
in_00020 <- intersect(hsa00020, matrix_genes)

cat("\n--- Overlap ---\n")
cat("hsa00010 genes in matrix:", length(in_00010), "/", length(hsa00010), "\n")
cat("hsa00030 genes in matrix:", length(in_00030), "/", length(hsa00030), "\n")
cat("hsa00020 genes in matrix:", length(in_00020), "/", length(hsa00020), "\n")

cat("\nMissing from hsa00010:", paste(setdiff(hsa00010, matrix_genes), collapse=", "), "\n")
cat("Missing from hsa00030:", paste(setdiff(hsa00030, matrix_genes), collapse=", "), "\n")
cat("Missing from hsa00020:", paste(setdiff(hsa00020, matrix_genes), collapse=", "), "\n")

cat("\nhsa00030 GENES (PPP) in matrix:\n")
cat(paste(in_00030, collapse=", "), "\n")

cat("\nhsa00020 GENES (TCA) in matrix:\n")
cat(paste(in_00020, collapse=", "), "\n")

matrix_only <- setdiff(matrix_genes, all_three)
cat("\nMatrix genes NOT in any of the three pathways (", length(matrix_only), "):\n")
cat(paste(matrix_only, collapse=", "), "\n")

# Save mapping
mapping <- data.frame(
  gene = all_three,
  in_hsa00010 = all_three %in% hsa00010,
  in_hsa00030 = all_three %in% hsa00030,
  in_hsa00020 = all_three %in% hsa00020,
  in_matrix = all_three %in% matrix_genes,
  stringsAsFactors = FALSE
)
rio::export(mapping, "results/tables/pathway_gene_membership.csv")
cat("\nPathway membership table saved.\n")
