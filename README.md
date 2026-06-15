# KIRP Glycolysis Transcriptomics

Análise transcriptômica da via glicólise/gliconeogênese em carcinoma renal papilar (KIRP), com foco em expressão diferencial, visualização por volcano plot e rede de interação proteína-proteína (PPI).

## Resumo científico

Este repositório organiza um fluxo reprodutível para investigar genes da via KEGG `hsa00010` (glicólise/gliconeogênese) em amostras de tecido renal normal e carcinoma renal papilar. O pipeline utiliza dados transcriptômicos TCGA/GTEx integrados via UCSC Xena, modelagem estatística com `limma`, correção de Benjamini-Hochberg, volcano plot e rede PPI com STRING.

## Justificativa

O carcinoma renal é uma neoplasia com forte componente de reprogramação metabólica. A glicólise/gliconeogênese é uma via relevante para investigar alterações transcriptômicas associadas à reorganização energética e biossintética tumoral. Este projeto prioriza transparência metodológica e rastreabilidade dos resultados, evitando análises não suportadas pelo fluxo principal.

## Objetivos

- Extrair genes humanos da via KEGG `hsa00010`.
- Preparar matriz TCGA/GTEx para comparação entre tecido renal normal e KIRP.
- Identificar genes diferencialmente expressos com `limma`.
- Gerar volcano plot em PNG com resolução mínima de 300 dpi.
- Construir rede PPI com STRING e exportar a figura final em PNG.

## Fonte dos dados

O pipeline espera uma matriz local `data/raw/kidney.tsv`, compatível com exportação UCSC Xena, contendo:

- colunas de metadados: `sample`, `primary_site`, `sample_type`, `study`, `TCGA_GTEX_main_category`;
- colunas de expressão gênica em escala `log2(norm_count + 1)`;
- categorias `GTEX Kidney` e `TCGA Kidney Papillary Cell Carcinoma`.

Por questões de tamanho, licença ou rastreabilidade, o arquivo bruto pode ser mantido fora do versionamento público quando necessário.

## Fluxo analítico

```text
KEGG hsa00010
  -> extração dos genes da glicólise/gliconeogênese
  -> importação TCGA/GTEx
  -> preparo da matriz de expressão
  -> limma: KIRP versus Normal
  -> Volcano Plot
  -> STRING
  -> Rede PPI
```

## Dependências

As principais dependências em R são:

- `KEGGREST`
- `dplyr`
- `stringr`
- `tibble`
- `rio`
- `limma`
- `ggplot2`
- `ggrepel`
- `STRINGdb`
- `igraph`
- `ggraph`

Arquivos de reprodutibilidade ficam em `environment/`.

## Estrutura do repositório

```text
KIRP_Glycolysis_Transcriptomics/
├── README.md
├── LICENSE
├── .gitignore
├── CITATION.cff
├── environment/
├── data/
├── scripts/
├── results/
├── docs/
└── output/
```

## Instruções de execução

1. Coloque a matriz `kidney.tsv` em `data/raw/`.
2. Abra R no diretório raiz do repositório.
3. Execute:

```r
source("scripts/run_pipeline.R")
```

## Exemplos de saída

As figuras finais são salvas exclusivamente em PNG:

- `results/figures/Volcano.png`
- `results/figures/PPI_network.png`

As tabelas principais são salvas em:

- `data/metadata/Hsa_genes.csv`
- `results/differential_expression/DEG_KIRP_vs_Normal.csv`
- `results/ppi/STRING_mapping.csv`
- `results/ppi/STRING_edges_high_confidence.csv`

## Escopo analítico

Este repositório mantém apenas:

- volcano plot;
- rede PPI.

Não são executadas análises de Kaplan-Meier, sobrevida, enriquecimento funcional, UpSet, dotplot, PPI 3D ou figuras suplementares.

## Contato

Ryan de Paulo Santos.

## Como citar

Use o arquivo `CITATION.cff` para citação automática pelo GitHub.
