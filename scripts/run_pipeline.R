required_packages <- c(
  "KEGGREST", "dplyr", "stringr", "tibble", "rio", "limma",
  "ggplot2", "ggrepel", "STRINGdb", "igraph", "ggraph", "umap"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Install missing packages before running the pipeline: ", paste(missing_packages, collapse = ", "))
}

scripts <- c(
  "01_download_data.R",
  "02_prepare_data.R",
  "02b_pca.R",
  "03_differential_expression.R",
  "04_volcano_plot.R",
  "05_ppi_network.R",
  "06_pathway_enrichment.R"
)

for (script in scripts) {
  script_path <- file.path("scripts", script)
  if (!file.exists(script_path)) {
    message("Skipping ", script, " (file not found)")
    next
  }
  message("\n>>> Running ", script)
  source(script_path, local = new.env(parent = globalenv()))
}

message("\nPipeline finished.")
