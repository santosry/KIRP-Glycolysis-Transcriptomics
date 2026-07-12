# run_all.R — Pipeline completo do KIRP transcriptomics
# Fluxo: kidney.tsv → DE → Volcano → PPI → Enrichment → Concordance
# Executar: Rscript run_all.R

message("═══════════════════════════════════════════")
message("  KIRP Central Carbon Metabolism Pipeline")
message("═══════════════════════════════════════════")

scripts <- c(
  "01_load_data.R",
  "02_differential_expression.R",
  "03_volcano_plots.R",
  "04_ppi_network.R",
  "05_functional_enrichment.R",
  "06_concordance_and_tables.R"
)

for(s in scripts) {
  message(sprintf("\n>>> Running: scripts/%s", s))
  t0 <- Sys.time()
  source(sprintf("scripts/%s", s))
  message(sprintf("<<< Done: %.1f min", difftime(Sys.time(), t0, units="mins")))
}

message("\n═══════════════════════════════════════════")
message("  PIPELINE COMPLETO")
message("═══════════════════════════════════════════")
message("Outputs: results/v3/")
message("  figures/ → PNG + 3D HTML")
message("  tables/  → DEGs, S1, camera, ORA, concordance")
message("  sessionInfo.txt")
