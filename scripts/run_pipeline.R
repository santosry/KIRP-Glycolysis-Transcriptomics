# run_pipeline.R
# 🍌 KIRP Glycolysis Transcriptomics — Pipeline completo
# Dois eixos independentes: (1) Transcriptoma global, (2) Via hsa00010

# ── Environment check ──
source("scripts/00_environment.R", local = new.env(parent = globalenv()))

# ── Pipeline scripts in order ──
pipeline_scripts <- c(
  "01_data_provenance.R",           # SHA256, scale check, gene IDs
  "02_prepare_data.R",              # Sample flow, metadata, expression matrix
  "03_sample_qc.R",                 # Boxplots, density, correlation
  "04_pca_umap.R",                  # PCA + UMAP + confounding audit
  "05_differential_expression_global.R",  # EIXO 1: DEG global + sensitivity
  "06_hsa00010_targeted_analysis.R",      # EIXO 2: KEGG hsa00010 audit + heatmap
  "07_ora_kegg.R",                  # ORA KEGG Up + Down
  "08_ora_reactome.R",              # ORA Reactome Up + Down
  "10_string_network.R",            # STRING + communities + centrality
  "11_integrative_analysis.R",      # Integration table
  "12_flowchart.R"                  # 🍌 Nano banana flowchart
)

for (script in pipeline_scripts) {
  script_path <- file.path("scripts", script)
  if (!file.exists(script_path)) {
    message("⚠ Skipping ", script, " (not found)")
    next
  }
  message("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  message("🍌 Running ", script)
  message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
  
  tryCatch({
    source(script_path, local = new.env(parent = globalenv()))
  }, error = function(e) {
    message("❌ ERROR in ", script, ": ", e$message)
    message("Continuing with next script...")
  })
}

message("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
message("🍌 PIPELINE COMPLETE! 🍌")
message("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
