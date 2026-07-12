setwd("C:/Users/oorie/OneDrive/Documentos/TRABALHOS/KIRP/KIRP-Glycolysis-Transcriptomics")
Sys.setenv(PATH = paste("C:\\Users\\oorie\\AppData\\Local\\Pandoc", Sys.getenv("PATH"), sep = ";"))
if(file.exists("manuscrito_kirp_tres_vias_FINAL_2026-07-11.Rmd")) {
  rmarkdown::render("manuscrito_kirp_tres_vias_FINAL_2026-07-11.Rmd", output_format = "pdf_document")
  cat("PDF rendered\n")
} else {
  cat("Rmd not found in repo. Copy from root first.\n")
}
