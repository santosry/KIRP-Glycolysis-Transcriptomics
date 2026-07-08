# run_pipeline.R — FAIL-FAST, com run_id, timestamp, manifest
# Interrompe imediatamente em qualquer erro. Não executa scripts dependentes.

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)

# ── Run ID ──
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
git_commit <- tryCatch(
  system("git rev-parse --short HEAD", intern = TRUE),
  error = function(e) "unknown"
)
run_id <- paste0("run_", timestamp, "_", git_commit)
message("Run ID: ", run_id)

# ── Create run directory ──
run_dir <- file.path(repo_root, "runs", run_id)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

# Copy all results to run directory at end
results_src <- file.path(repo_root, "results")
data_src    <- file.path(repo_root, "data")

# ── Manifest ──
manifest <- data.frame(
  run_id        = character(),
  git_commit    = character(),
  script        = character(),
  start_time    = character(),
  end_time      = character(),
  status        = character(),
  output_files  = character(),
  error_message = character(),
  stringsAsFactors = FALSE
)

write_manifest <- function() {
  write.csv(manifest, file.path(run_dir, "pipeline_manifest.csv"), row.names = FALSE)
}

# ── Required environment check ──
message("\n========== 00_environment.R ==========")
start_t <- Sys.time()
tryCatch({
  source("scripts/00_environment.R", local = new.env(parent = globalenv()))
  status <- "SUCCESS"; err <- ""
}, error = function(e) {
  message("FATAL: Environment check failed.")
  message(e$message)
  status <- "FAILED"; err <- e$message
  manifest <<- rbind(manifest, data.frame(
    run_id = run_id, git_commit = git_commit,
    script = "00_environment.R",
    start_time = as.character(start_t), end_time = as.character(Sys.time()),
    status = status, output_files = "", error_message = err,
    stringsAsFactors = FALSE
  ))
  write_manifest()
  stop("Pipeline aborted: environment check failed. Fix missing packages and retry.")
})
manifest <- rbind(manifest, data.frame(
  run_id = run_id, git_commit = git_commit,
  script = "00_environment.R",
  start_time = as.character(start_t), end_time = as.character(Sys.time()),
  status = status, output_files = "", error_message = err,
  stringsAsFactors = FALSE
))

# ── Optional: download data ──
download_script <- "scripts/00_download_data.R"
if (file.exists(download_script)) {
  message("\n========== 00_download_data.R ==========")
  start_t <- Sys.time()
  tryCatch({
    source(download_script, local = new.env(parent = globalenv()))
    status <- "SUCCESS"; err <- ""
  }, error = function(e) {
    message("WARNING: Download script failed (may need manual download): ", e$message)
    status <- "WARNING"; err <- e$message
  })
  manifest <- rbind(manifest, data.frame(
    run_id = run_id, git_commit = git_commit,
    script = "00_download_data.R",
    start_time = as.character(start_t), end_time = as.character(Sys.time()),
    status = status, output_files = "", error_message = err,
    stringsAsFactors = FALSE
  ))
}

# ── Pipeline scripts in strict order ──
pipeline_scripts <- c(
  "01_data_provenance.R",
  "02_prepare_data.R",
  "03_sample_qc.R",
  "04_pca_umap.R",
  "05_differential_expression_global.R",
  "06_hsa00010_targeted_analysis.R",
  "07_ora_kegg.R",
  "08_ora_reactome.R",
  "09_rank_based_enrichment.R",
  "10_string_network.R",
  "11_integrative_analysis.R",
  "12_flowchart.R",
  "13_generate_manuscript_results.R"
)

all_success <- TRUE
for (script in pipeline_scripts) {
  script_path <- file.path("scripts", script)
  cat("\n========================================\n")
  cat("RUNNING:", script, "\n")
  cat("========================================\n")
  start_t <- Sys.time()
  
  if (!file.exists(script_path)) {
    message("SKIPPING ", script, " (file not found)")
    manifest <- rbind(manifest, data.frame(
      run_id = run_id, git_commit = git_commit,
      script = script, start_time = as.character(start_t),
      end_time = as.character(Sys.time()),
      status = "SKIPPED", output_files = "", error_message = "File not found",
      stringsAsFactors = FALSE
    ))
    next
  }
  
  result <- tryCatch({
    source(script_path, local = new.env(parent = globalenv()))
    "SUCCESS"
  }, error = function(e) {
    message("\n========================================")
    message("FATAL ERROR in ", script)
    message("========================================")
    message(e$message)
    return(paste("FAILED:", e$message))
  })
  
  end_t <- Sys.time()
  status <- if (result == "SUCCESS") "SUCCESS" else "FAILED"
  
  manifest <- rbind(manifest, data.frame(
    run_id = run_id, git_commit = git_commit,
    script = script,
    start_time = as.character(start_t), end_time = as.character(end_t),
    status = status, output_files = "", error_message = if(result != "SUCCESS") result else "",
    stringsAsFactors = FALSE
  ))
  
  if (status == "FAILED") {
    all_success <- FALSE
    message("\nPipeline ABORTED after failure in ", script)
    break
  }
}

# ── Finalize ──
session_file <- file.path(run_dir, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), session_file)

manifest$output_files <- sapply(pipeline_scripts, function(s) {
  # Collect output files that exist
  script_num <- sub("_.*", "", s)
  pattern <- paste0("results/.*", script_num, ".*")
  files <- list.files(file.path(repo_root, "results"), recursive = TRUE, full.names = TRUE)
  paste(files, collapse = "; ")
})

write_manifest()

# ── Copy results to run directory ──
if (dir.exists(results_src)) {
  file.copy(results_src, run_dir, recursive = TRUE, overwrite = TRUE)
}

# ── Final status ──
message("\n========================================")
if (all_success) {
  message("PIPELINE COMPLETE — ALL STEPS SUCCESSFUL")
  message("Run ID: ", run_id)
  message("Manifest: runs/", run_id, "/pipeline_manifest.csv")
} else {
  message("PIPELINE INCOMPLETE — errors occurred")
  message("Run ID: ", run_id)
  message("Check runs/", run_id, "/pipeline_manifest.csv for details")
}
message("========================================")
