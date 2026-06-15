required_packages <- c(
  "KEGGREST", "dplyr", "stringr", "tibble", "rio", "limma",
  "ggplot2", "ggrepel", "STRINGdb", "igraph", "ggraph"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Install missing packages before running the pipeline: ", paste(missing_packages, collapse = ", "))
}

scripts <- c(
  "01_download_data.R",
  "02_prepare_data.R",
  "03_differential_expression.R",
  "04_volcano_plot.R",
  "05_ppi_network.R"
)

for (script in scripts) {
  message("\n>>> Running ", script)
  source(file.path("scripts", script), local = new.env(parent = globalenv()))
}

message("\nPipeline finished.")
