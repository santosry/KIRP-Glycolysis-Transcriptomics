# Workflow

```text
KEGG
  -> TCGA/GTEx
  -> limma
  -> Volcano Plot
  -> STRING
  -> PPI
```

## 1. KEGG

O script `01_download_data.R` consulta a via `hsa00010` no KEGG usando `KEGGREST::keggGet()` e extrai os símbolos gênicos humanos relacionados à glicólise/gliconeogênese.

## 2. TCGA/GTEx

O script `02_prepare_data.R` importa `data/raw/kidney.tsv`, uma matriz TCGA/GTEx em escala `log2(norm_count + 1)`. As amostras são classificadas como:

- `Normal`: `GTEX Kidney`;
- `KIRP`: `TCGA Kidney Papillary Cell Carcinoma`.

## 3. limma

O script `03_differential_expression.R` ajusta um modelo linear com `limma`, define o contraste `KIRP - Normal`, aplica `eBayes()` e corrige os valores de p pelo método de Benjamini-Hochberg.

## 4. Volcano Plot

O script `04_volcano_plot.R` gera `results/figures/Volcano.png`, usando FDR < 0,05 e `|logFC| > 1` como critérios visuais e analíticos.

## 5. STRING

O script `05_ppi_network.R` mapeia os DEGs para identificadores STRING, recupera interações proteína-proteína e mantém apenas arestas com `combined_score >= 700`.

## 6. PPI

A rede é convertida para `igraph`, o maior componente conexo é visualizado com `ggraph` e a figura final é salva como `results/figures/PPI_network.png`.
