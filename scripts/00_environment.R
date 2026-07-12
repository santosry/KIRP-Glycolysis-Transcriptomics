# 00_environment.R
# Verifica e registra ambiente. stop() se pacotes ausentes.
# v3.1.0 — apenas pacotes essenciais para o pipeline v3

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
dir.create(file.path(repo_root, "environment"), recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "dplyr", "tibble", "stringr", "rio", "digest",
  "limma",
  "ggplot2", "ggrepel", "pheatmap", "RColorBrewer",
  "KEGGREST",
  "clusterProfiler", "org.Hs.eg.db",
  "igraph",
  "plotly", "htmlwidgets",
  "tidyr"
)

message("Checking required packages...")
missing <- c()
for (pkg in required_packages) {
  has_it <- requireNamespace(pkg, quietly = TRUE)
  if (!has_it) missing <- c(missing, pkg)
}

if (length(missing) > 0) {
  msg <- paste0(
    "FATAL: Missing required packages: ", paste(missing, collapse = ", "), "\n",
    "Install with:\n",
    "  install.packages(c('", paste(setdiff(missing, c("limma","clusterProfiler","org.Hs.eg.db","KEGGREST")), collapse = "','"), "'))\n",
    "  BiocManager::install(c('", paste(intersect(missing, c("limma","clusterProfiler","org.Hs.eg.db","KEGGREST")), collapse = "','"), "'))"
  )
  stop(msg)
}

message("All ", length(required_packages), " required packages available.")
message("R version: ", R.version.string)

if (requireNamespace("BiocManager", quietly = TRUE)) {
  message("Bioconductor version: ", as.character(BiocManager::version()))
}

# Save package versions
pkg_versions <- data.frame(
  package = required_packages,
  version = sapply(required_packages, function(p) as.character(packageVersion(p))),
  stringsAsFactors = FALSE
)
rio::export(pkg_versions, file.path(repo_root, "environment", "packages.csv"))
message("Environment recorded.")
