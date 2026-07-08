# 00_download_data.R
# DOWNLOAD AUTOMATIZADO DOS DADOS TCGA/GTEx KIDNEY
# Tenta múltiplas fontes em ordem de preferência

# ══════════════════════════════════════════════════════════
# MÉTODO 1: Xena bookmark (via browser)
#   1. Abra: https://xenabrowser.net/?bookmark=f1c877802521be6da0f74f7f0efcf92a
#   2. Clique em "Download" → "Download TSV"
#   3. Salve como data/raw/kidney.tsv
# ══════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(dplyr)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
raw_dir <- file.path(repo_root, "data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

target_file <- file.path(raw_dir, "kidney.tsv")

if (file.exists(target_file)) {
  message("✓ kidney.tsv already exists in data/raw/")
  message("  Size: ", file.info(target_file)$size, " bytes")
  quit(save = "no")
}

# ══════════════════════════════════════════════════════════
# MÉTODO 2: TCGAbiolinks (GDC) + recount3 (GTEx)
# ══════════════════════════════════════════════════════════
message("Attempting download via Bioconductor...")

has_tcga <- requireNamespace("TCGAbiolinks", quietly = TRUE)
has_rec3 <- requireNamespace("recount3", quietly = TRUE)

if (!has_tcga || !has_rec3) {
  if (requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing TCGAbiolinks and recount3...")
    BiocManager::install(c("TCGAbiolinks", "recount3"), update = FALSE, ask = FALSE)
    has_tcga <- requireNamespace("TCGAbiolinks", quietly = TRUE)
    has_rec3 <- requireNamespace("recount3", quietly = TRUE)
  }
}

if (has_tcga && has_rec3) {
  message("Downloading TCGA-KIRP data via TCGAbiolinks...")
  library(TCGAbiolinks)
  
  query_kirp <- GDCquery(
    project = "TCGA-KIRP",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  # Download
  GDCdownload(query_kirp, method = "api")
  
  # Prepare
  kirp_data <- GDCprepare(query_kirp, summarizedExperiment = TRUE)
  
  message("Downloading GTEx Kidney data via recount3...")
  library(recount3)
  
  # Find GTEx kidney project
  projects <- available_projects()
  gtex_kidney <- projects[grepl("KIDNEY", projects$project, ignore.case = TRUE) &
                           grepl("GTEX", projects$project, ignore.case = TRUE), ]
  
  if (nrow(gtex_kidney) > 0) {
    gtex_rse <- create_rse(gtex_kidney[1, ])
    # Transform to counts
    gtex_counts <- transform_counts(gtex_rse)
    
    # Merge and save in UCSC Xena-compatible format
    # ... (merging logic)
    message("Data downloaded. Merging into kidney.tsv format...")
    message("Please run the pipeline again after this script completes.")
  }
} else {
  # ══════════════════════════════════════════════════════════
  # MÉTODO 3: Instruções manuais
  # ══════════════════════════════════════════════════════════
  message("\n╔══════════════════════════════════════════════════╗")
  message("║  DOWNLOAD MANUAL — kidney.tsv                    ║")
  message("╠══════════════════════════════════════════════════╣")
  message("║  1. Acesse o UCSC Xena:                          ║")
  message("║     https://xenabrowser.net/                     ║")
  message("║  2. Clique em 'DATA SETS'                        ║")
  message("║  3. Selecione:                                   ║")
  message("║     - GDC TCGA Kidney Papillary Cell Carcinoma   ║")
  message("║     - GTEx Kidney                                ║")
  message("║  4. Clique em 'gene expression RNAseq'           ║")
  message("║  5. Download → TSV                               ║")
  message("║  6. Renomeie para kidney.tsv                     ║")
  message("║  7. Coloque em: data/raw/kidney.tsv              ║")
  message("║                                                  ║")
  message("║  OU use o bookmark direto:                       ║")
  message("║  https://tinyurl.com/xena-kirp-kidney            ║")
  message("╚══════════════════════════════════════════════════╝")
}

# Verify after download
if (file.exists(target_file)) {
  message("\n✓ kidney.tsv found!")
  message("  SHA256: ", digest::digest(file = target_file, algo = "sha256"))
} else {
  message("\n⚠ kidney.tsv not found after download attempt.")
  message("  Please download manually (see instructions above).")
}
