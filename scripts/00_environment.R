suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(rio)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)

# ── Required packages for full pipeline ──
required_packages <- c(
  # Data & QC
  "KEGGREST", "dplyr", "stringr", "tibble", "rio",
  # Differential expression
  "limma", "edgeR",
  # Visualization
  "ggplot2", "ggrepel", "pheatmap", "RColorBrewer",
  # Dimensionality reduction
  "umap",
  # Enrichment
  "clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot",
  # Network
  "STRINGdb", "igraph", "ggraph",
  # Misc
  "digest", "tools"
)

message("Checking required packages...")
missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  message("Missing packages: ", paste(missing, collapse = ", "))
  message("Install with:\n  install.packages(c('", paste(setdiff(missing, c("limma","edgeR","clusterProfiler","org.Hs.eg.db","ReactomePA","enrichplot","STRINGdb")), collapse = "','"), "'))\n  BiocManager::install(c('", paste(intersect(missing, c("limma","edgeR","clusterProfiler","org.Hs.eg.db","ReactomePA","enrichplot","STRINGdb")), collapse = "','"), "'))")
} else {
  message("All packages available.")
}

# ── Environment info ──
dir.create(file.path(repo_root, "environment"), recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(repo_root, "environment", "sessionInfo.txt"))

pkg_info <- as.data.frame(installed.packages()[, c("Package","Version")])
pkg_info <- pkg_info[pkg_info$Package %in% required_packages, ]
rio::export(pkg_info, file.path(repo_root, "environment", "packages.csv"))

message("Environment recorded.")
message("R version: ", R.version.string)
message("Bioconductor: ", as.character(BiocManager::version()))
