# KIRP Glycolysis Transcriptomics

Análise transcriptômica da via glicólise/gliconeogênese em carcinoma renal papilar (KIRP), com subsídios moleculares para a prática do enfermeiro em genética, genômica e saúde de precisão.

## Resumo científico

Este repositório organiza um fluxo reprodutível para investigar genes da via KEGG `hsa00010` (glicólise/gliconeogênese) em amostras de tecido renal normal e carcinoma renal papilar. O *pipeline* utiliza dados transcriptômicos TCGA/GTEx integrados via UCSC Xena, modelagem estatística com *limma*, correção de Benjamini-Hochberg, *volcano plot* e rede PPI com STRING. O manuscrito derivado discute como o enfermeiro que atua em genética, genômica e saúde de precisão pode interpretar esses achados moleculares para qualificar sua prática clínica.

## Autores

- **Kamila da Conceição Loureiro** — Faculdade de Medicina de Campos, Campos dos Goytacazes, RJ, Brasil.
- **Ryan de Paulo Santos** — Instituto Federal de Educação, Ciência e Tecnologia Fluminense Campus Campos Guarus, Campos dos Goytacazes, RJ, Brasil.
- **Letícia Maria Dias Freitas** — Escola Técnica Estadual João Barcelos Martins, Campos dos Goytacazes, RJ, Brasil.

## Manuscrito

O manuscrito completo está disponível em `manuscript/manuscrito_kirp_glicolise.Rmd` e `manuscript/manuscrito_kirp_glicolise.pdf`.

**Título:** Transcriptômica da glicólise no carcinoma renal papilar: subsídios moleculares para a prática do enfermeiro em genética, genômica e saúde de precisão

**Palavras-chave:** Carcinoma de Células Renais; Glicólise; Genômica; Enfermagem de Precisão; Enfermagem Genética.

## Justificativa

O carcinoma renal é uma neoplasia com forte componente de reprogramação metabólica. A glicólise/gliconeogênese é uma via relevante para investigar alterações transcriptômicas associadas à reorganização energética e biossintética tumoral. Paralelamente, a enfermagem contemporânea testemunha a consolidação da enfermagem em genética e genômica e da enfermagem de precisão, campos que demandam alfabetização molecular para interpretação de evidências ômicas. Este projeto prioriza transparência metodológica, rastreabilidade dos resultados e tradução do conhecimento para a prática clínica da enfermagem.

## Objetivos

- Extrair genes humanos da via KEGG `hsa00010`.
- Preparar matriz TCGA/GTEx para comparação entre tecido renal normal e KIRP.
- Identificar genes diferencialmente expressos com *limma*.
- Gerar *volcano plot* em PNG com resolução mínima de 300 dpi.
- Construir rede PPI com STRING e exportar a figura final em PNG.
- Discutir as implicações translacionais dos achados para a enfermagem em genética, genômica e saúde de precisão.

## Fonte dos dados

O *pipeline* espera uma matriz local `data/raw/kidney.tsv`, compatível com exportação UCSC Xena, contendo:

- colunas de metadados: `sample`, `primary_site`, `sample_type`, `study`, `TCGA_GTEX_main_category`;
- colunas de expressão gênica em escala *log2(norm\_count + 1)*;
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
KIRP-Glycolysis-Transcriptomics/
├── README.md
├── LICENSE
├── .gitignore
├── CITATION.cff
├── environment/
├── data/
├── scripts/
├── results/
├── docs/
├── manuscript/
├── output/
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

- *volcano plot*;
- rede PPI.

Não são executadas análises de Kaplan-Meier, sobrevida, enriquecimento funcional, UpSet, *dotplot*, PPI 3D ou figuras suplementares.

## Declaração de uso de Inteligência Artificial

Em conformidade com a Portaria CNPq nº 2.664, de 22 de dezembro de 2026, declara-se que ferramentas de inteligência artificial generativa foram utilizadas como recursos de apoio durante a elaboração do manuscrito, da documentação e da organização deste repositório. Especificamente, foram empregados ChatGPT-5.5 (OpenAI), Codex (OpenAI) e DeepSeek-v4-pro (Hangzhou DeepSeek Artificial Intelligence) nas seguintes etapas: (i) revisão gramatical, estilística e de coesão textual; (ii) sugestão de estrutura argumentativa para a seção de discussão do manuscrito; (iii) verificação de consistência entre as referências citadas e o conteúdo discutido; (iv) formatação de referências conforme o estilo Vancouver; e (v) organização da documentação do repositório. Em nenhum momento, as ferramentas de IA substituíram o julgamento científico dos autores, a análise crítica dos resultados, a interpretação dos dados transcriptômicos ou a responsabilidade pelo conteúdo final. Todos os autores revisaram, editaram e aprovam integralmente o conteúdo, assumindo total responsabilidade intelectual pelo repositório e pelo manuscrito.

## Contato

ryan.paulo@gsuite.iff.edu.br

## Como citar

Use o arquivo `CITATION.cff` para citação automática pelo GitHub.
