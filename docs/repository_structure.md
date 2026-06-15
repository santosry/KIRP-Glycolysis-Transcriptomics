# Repository Structure

```text
KIRP_Glycolysis_Transcriptomics/
├── README.md
├── LICENSE
├── .gitignore
├── CITATION.cff
├── environment/
│   ├── sessionInfo.txt
│   └── packages.csv
├── data/
│   ├── raw/
│   ├── processed/
│   └── metadata/
├── scripts/
│   ├── 01_download_data.R
│   ├── 02_prepare_data.R
│   ├── 03_differential_expression.R
│   ├── 04_volcano_plot.R
│   ├── 05_ppi_network.R
│   └── run_pipeline.R
├── results/
│   ├── differential_expression/
│   ├── figures/
│   ├── ppi/
│   └── tables/
├── docs/
│   ├── workflow.md
│   ├── methods_summary.md
│   └── repository_structure.md
└── output/
```

## Versionamento

O manuscrito local em R Markdown e seu PDF derivado não devem ser versionados. O `.gitignore` bloqueia `*.Rmd`, `*.pdf`, `output/`, `manuscript/`, caches, logs e arquivos temporários.

## Dados brutos

O diretório `data/raw/` é destinado ao arquivo local `kidney.tsv`. Dependendo das permissões de redistribuição e tamanho, esse arquivo pode permanecer fora do GitHub.
