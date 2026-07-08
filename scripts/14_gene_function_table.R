# 14_gene_function_table.R
# Consulta STRING para anotações funcionais dos 31 DEGs

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(rio); library(STRINGdb)
})

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
deg_file <- file.path(repo_root, "results", "differential_expression", "DEG_global.csv")

if (!file.exists(deg_file)) stop("Run 05 first.")

deg <- rio::import(deg_file) |> filter(regulation %in% c("Up", "Down"))

# STRING functional annotations
options(timeout = 5000000)
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 0, input_directory = "")

# Map all DEGs to STRING
mapped <- string_db$map(as.data.frame(deg |> select(gene_id) |> rename(gene_symbol = gene_id)), 
                         "gene_symbol", removeUnmappedRows = TRUE)
mapped <- mapped |> filter(!is.na(STRING_id))

# Get functional annotations for each protein
func_table <- tibble(
  gene_symbol = character(),
  preferred_name = character(),
  protein_size = integer(),
  annotation = character(),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(mapped)) {
  sid <- mapped$STRING_id[i]
  gs <- mapped$gene_symbol[i]
  
  # Try to get annotation
  annot <- tryCatch({
    string_db$get_annotations(sid)
  }, error = function(e) NULL)
  
  if (!is.null(annot) && nrow(annot) > 0) {
    func_table <- rbind(func_table, tibble(
      gene_symbol = gs,
      preferred_name = annot$preferred_name[1],
      protein_size = annot$protein_size[1],
      annotation = annot$annotation[1]
    ))
  } else {
    func_table <- rbind(func_table, tibble(
      gene_symbol = gs,
      preferred_name = NA_character_,
      protein_size = NA_integer_,
      annotation = "Not available"
    ))
  }
}

# Merge with DE results
func_table <- func_table |> 
  left_join(deg |> select(gene_id, logFC, adj.P.Val, regulation), by = c("gene_symbol" = "gene_id")) |>
  arrange(desc(abs(logFC)))

rio::export(func_table, file.path(repo_root, "results", "tables", "gene_function_table.csv"))

# Print summary
message("Gene function table: ", nrow(func_table), " genes annotated")
for (i in 1:min(10, nrow(func_table))) {
  r <- func_table[i, ]
  message(sprintf("  %-10s [%s] %s", r$gene_symbol, r$regulation, substr(r$annotation, 1, 80)))
}

message("\nSaved results/tables/gene_function_table.csv")
