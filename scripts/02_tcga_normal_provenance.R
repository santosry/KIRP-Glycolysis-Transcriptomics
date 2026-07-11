# 02_tcga_normal_provenance.R
# Characterize all 129 TCGA_Normal samples
# Determine which TCGA projects they belong to

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(rio); library(stringr)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
raw_file <- file.path(repo_root, "data", "raw", "kidney.tsv")

dir.create(file.path(repo_root, "results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "results", "sensitivity"), recursive = TRUE, showWarnings = FALSE)

kidney <- rio::import(raw_file)

# Filter TCGA Normal samples (blank category, TCGA study, Solid Tissue Normal)
tcga_norm <- kidney |>
  filter(`_study` == "TCGA",
         `_sample_type` == "Solid Tissue Normal",
         is.na(TCGA_GTEX_main_category) | TCGA_GTEX_main_category == "")

message("TCGA Normal samples (Solid Tissue Normal, blank category): ", nrow(tcga_norm))

# Parse TCGA barcodes
# Format: TCGA-XX-YYYY-ZZ where:
#   XX = tissue source site
#   YYYY = participant
#   ZZ = sample type (01=primary tumor, 11=solid tissue normal, etc.)
# The TCGA project (KIRP, KIRC, KICH) is NOT encoded in the barcode itself
# We need to query the GDC API or use known project mappings

tcga_norm <- tcga_norm |>
  mutate(
    sample_id = sample,
    barcode_parts = str_split(sample, "-"),
    tss = sapply(barcode_parts, function(x) if(length(x) >= 2) x[2] else NA),
    participant = sapply(barcode_parts, function(x) if(length(x) >= 3) x[3] else NA),
    sample_code = sapply(barcode_parts, function(x) if(length(x) >= 4) {
      sc <- substr(x[4], 1, 2); sc
    } else NA),
    vial = sapply(barcode_parts, function(x) if(length(x) >= 4) x[4] else NA)
  )

# Sample type codes: 01=Primary Solid Tumor, 11=Solid Tissue Normal
message("\nSample type codes in TCGA Normals:")
print(table(tcga_norm$sample_code))

# Try to infer project from the sample ID prefix
# TCGA-KIRP samples typically have barcodes starting with TCGA-XX
# where XX maps to specific projects
# Without GDC API, we can check if the participant ID overlaps with KIRP tumor participants

kirp_tumors <- kidney |>
  filter(TCGA_GTEX_main_category == "TCGA Kidney Papillary Cell Carcinoma")

kirp_participants <- kirp_tumors |>
  mutate(participant = sapply(str_split(sample, "-"), function(x) if(length(x) >= 3) x[3] else NA)) |>
  pull(participant) |> unique()

message("\nKIRP tumor participants: ", length(kirp_participants))

# Check overlap: how many TCGA Normals share participants with KIRP tumors?
tcga_norm$matched_to_kirp <- tcga_norm$participant %in% kirp_participants
n_kirp_matched <- sum(tcga_norm$matched_to_kirp)
n_other <- nrow(tcga_norm) - n_kirp_matched

message("TCGA Normals matched to KIRP participants: ", n_kirp_matched)
message("TCGA Normals NOT matched to KIRP participants: ", n_other, 
        " (likely from KIRC, KICH, or other projects)")

# Build provenance table
provenance_table <- tcga_norm |>
  select(sample_id = sample, participant, sample_code, vial, matched_to_kirp) |>
  mutate(
    project_inferred = ifelse(matched_to_kirp, "TCGA-KIRP (likely)", "TCGA-OTHER (KIRC/KICH)"),
    sample_type = "Solid Tissue Normal",
    tissue = "Kidney",
    comparator_suitability = ifelse(matched_to_kirp, 
                                     "Suitable: KIRP-matched normal", 
                                     "Caution: may not be KIRP-specific normal")
  )

rio::export(provenance_table, file.path(repo_root, "results", "sensitivity", "tcga_normal_provenance.csv"))

# Summary
message("\n=== TCGA NORMAL PROVENANCE ===")
message("Total: ", nrow(tcga_norm))
message("Matched to KIRP participants: ", n_kirp_matched)
message("Not matched (other renal projects): ", n_other)
message("Sample codes: all ", paste(unique(tcga_norm$sample_code), collapse=", "))
message("\nImplications:")
message("  - ", n_kirp_matched, " TCGA_Normal samples share participant IDs with KIRP tumors")
message("    These may represent paired tumor-normal samples from KIRP patients")
message("  - ", n_other, " samples do not match KIRP participants")
message("    These likely originate from KIRC, KICH, or other TCGA kidney projects")

# Check for actual paired samples
if (n_kirp_matched > 0) {
  # Find participants with both tumor and normal
  kirp_tumor_participants <- kirp_tumors |>
    mutate(participant = sapply(str_split(sample, "-"), function(x) if(length(x) >= 3) x[3] else NA))
  
  normal_with_pair <- tcga_norm |>
    filter(matched_to_kirp) |>
    select(sample, participant)
  
  # Which participants have both?
  both <- intersect(
    unique(normal_with_pair$participant),
    unique(kirp_tumor_participants$participant)
  )
  message("\nParticipants with BOTH tumor and normal samples: ", length(both))
  
  # Export paired sample list
  pairs <- tibble(
    participant = both,
    has_tumor = TRUE,
    has_normal = TRUE
  )
  rio::export(pairs, file.path(repo_root, "results", "sensitivity", "paired_samples.csv"))
  
  message("  → ", length(both), " potential paired tumor-normal samples available for paired analysis")
}

message("\n✓ TCGA Normal provenance complete.")
