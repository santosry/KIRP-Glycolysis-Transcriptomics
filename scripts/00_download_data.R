# 00_download_data.R
# INSTRUÇÕES PARA DOWNLOAD MANUAL DO kidney.tsv
# O arquivo NÃO pode ser baixado programaticamente.
# Deve ser obtido manualmente via UCSC Xena.

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(repo_root, "data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
target_file <- file.path(raw_dir, "kidney.tsv")

if (file.exists(target_file) && file.info(target_file)$size > 1000000) {
  message("✓ kidney.tsv found (", round(file.info(target_file)$size/1e6, 1), " MB)")
  message("  SHA256: ", digest::digest(file = target_file, algo = "sha256"))
  quit(save = "no")
}

cat("
╔══════════════════════════════════════════════════════════════╗
║         DOWNLOAD MANUAL DO kidney.tsv                       ║
╠══════════════════════════════════════════════════════════════╣
║                                                            ║
║  1. Abra o bookmark no navegador:                          ║
║     https://xenabrowser.net/?bookmark=f1c877802521be6da0f74f7f0efcf92a
║                                                            ║
║  2. O Xena carregará o dataset combinado:                  ║
║     - GDC TCGA Kidney Papillary Cell Carcinoma (KIRP)      ║
║     - GTEx Kidney                                          ║
║                                                            ║
║  3. Clique em 'Download' → 'Download TSV'                  ║
║                                                            ║
║  4. Renomeie o arquivo para: kidney.tsv                    ║
║                                                            ║
║  5. Coloque o arquivo em:                                  ║
║     data/raw/kidney.tsv                                    ║
║                                                            ║
║  6. Execute novamente: source('scripts/run_pipeline.R')    ║
║                                                            ║
╚══════════════════════════════════════════════════════════════╝
")

stop("kidney.tsv not found. Please download manually from UCSC Xena bookmark (see instructions above).")
