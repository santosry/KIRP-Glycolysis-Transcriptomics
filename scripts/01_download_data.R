suppressPackageStartupMessages({
  library(KEGGREST)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)

dir.create(file.path(repo_root, "data", "metadata"), recursive = TRUE, showWarnings = FALSE)

pathway_id <- "hsa00010"
pw <- keggGet(pathway_id)
stopifnot(length(pw) == 1L)

genes_raw <- pw[[1]]$GENE
if (is.null(genes_raw) || length(genes_raw) == 0L) {
  stop("No GENE field returned by KEGG for ", pathway_id)
}

gene_symbols <- genes_raw[seq(2, length(genes_raw), by = 2)] |>
  str_remove("\\s*\\[.*$") |>
  str_trim() |>
  str_extract("^[^;]+") |>
  str_trim() |>
  unique() |>
  sort()

genes_via <- tibble(
  pathway_id = pathway_id,
  gene_symbol = gene_symbols
)

rio::export(genes_via, file.path(repo_root, "data", "metadata", "Hsa_genes.csv"))
message("Saved KEGG gene list: data/metadata/Hsa_genes.csv")
