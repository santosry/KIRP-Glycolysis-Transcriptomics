# KIRP Glycolysis Transcriptomics

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R 4.6.0](https://img.shields.io/badge/R-4.6.0-276DC3)](https://www.r-project.org/)

**Análise transcriptômica reprodutível da via glicólise/gliconeogênese em carcinoma renal papilar (KIRP)**

---

## Sumário

- [Sobre o projeto](#sobre-o-projeto)
- [Autores](#autores)
- [Fundamentação científica](#fundamentação-científica)
- [Desenho do estudo](#desenho-do-estudo)
- [Fluxo analítico](#fluxo-analítico)
- [Estrutura do repositório](#estrutura-do-repositório)
- [Dependências](#dependências)
- [Fonte dos dados](#fonte-dos-dados)
- [Como executar o pipeline](#como-executar-o-pipeline)
- [Descrição detalhada dos scripts](#descrição-detalhada-dos-scripts)
- [Saídas geradas](#saídas-geradas)
- [Escopo e limitações analíticas](#escopo-e-limitações-analíticas)
- [Reprodutibilidade](#reprodutibilidade)
- [Declaração de uso de Inteligência Artificial](#declaração-de-uso-de-inteligência-artificial)
- [Licença](#licença)
- [Contato](#contato)
- [Como citar](#como-citar)

---

## Sobre o projeto

Este repositório organiza um *pipeline* computacional completo e reprodutível para investigar a expressão de genes da via glicólise/gliconeogênese (KEGG `hsa00010`) em amostras de carcinoma renal papilar (*Kidney Renal Papillary Cell Carcinoma* — KIRP) e tecido renal normal.

O projeto integra:

- **Dados transcriptômicos públicos** do TCGA e GTEx, acessados via UCSC Xena;
- **Análise de expressão diferencial** com modelagem estatística robusta (*limma*);
- **Visualização** por *volcano plot* em alta resolução;
- **Rede de interação proteína-proteína (PPI)** construída com STRING;
- **Discussão translacional** voltada à alfabetização genômica do enfermeiro em saúde de precisão.

Cada etapa é executada por um script independente, permitindo auditoria, reexecução seletiva e adaptação a outras vias metabólicas ou tipos tumorais.

---

## Autores

| Autor | Afiliação |
|-------|-----------|
| **Kamila da Conceição Loureiro** | Faculdade de Medicina de Campos, Campos dos Goytacazes, RJ, Brasil |
| **Ryan de Paulo Santos** | Instituto Federal de Educação, Ciência e Tecnologia Fluminense *Campus* Campos Guarus, Campos dos Goytacazes, RJ, Brasil |
| **Letícia Maria Dias Freitas** | Escola Técnica Estadual João Barcelos Martins, Campos dos Goytacazes, RJ, Brasil |
| **Ivine** | Universidade Federal do Rio de Janeiro, *Campus* Macaé, Macaé, RJ, Brasil |
| **Maria Eduarda Pecly** | Instituto Federal de Educação, Ciência e Tecnologia Fluminense *Campus* Campos Guarus, Campos dos Goytacazes, RJ, Brasil |

---

## Fundamentação científica

### Contexto biológico

O carcinoma renal papilar é o segundo subtipo histológico mais frequente de carcinoma de células renais (10–15% dos casos). Diferentemente do carcinoma de células claras (KIRC), o KIRP exibe paisagem genética distinta — mutações em *MET*, ganhos dos cromossomos 7 e 17, perda do cromossomo Y — e reprogramação metabólica que envolve as vias glicolítica, lipídica e do ciclo do ácido tricarboxílico.

A glicólise ocupa posição central nessa reprogramação. O **efeito Warburg** — manutenção de elevada taxa glicolítica mesmo em presença de oxigênio — sustenta a proliferação celular tumoral por meio do fornecimento de ATP, intermediários biossintéticos (pentoses-fosfato, aminoácidos, lipídios) e equivalentes redox (NADH, NADPH).

### Por que a via `hsa00010`?

A via KEGG `hsa00010` (glicólise/gliconeogênese) foi escolhida como objeto de estudo porque:

1. **É uma via metabólica central**, evolutivamente conservada e bem anotada;
2. **Contém 65 genes humanos** com funções enzimáticas e regulatórias bem caracterizadas;
3. **Permite leitura em nível de sistemas**: os genes não são independentes, mas organizam-se em módulos funcionais (fase preparatória, fase de pagamento, gliconeogênese, metabolismo de aldeídos);
4. **Conecta-se ao efeito Warburg**: a reprogramação glicolítica é uma *hallmark* do câncer. No efeito Warburg, células tumorais mantêm elevada taxa glicolítica mesmo em presença de oxigênio (glicólise aeróbica), desviando o piruvato para lactato em vez de oxidá-lo no ciclo do ácido tricarboxílico. Esse desvio metabólico — aparentemente paradoxal, dado o menor rendimento de ATP por glicose (2 ATP *vs.* ~36 ATP na fosforilação oxidativa) — é compensado pela alta velocidade da glicólise e, sobretudo, pelo fornecimento de intermediários biossintéticos essenciais à proliferação celular: ribose-5-fosfato (via pentose-fosfato), serina e glicina (via 3-fosfoglicerato), e equivalentes redox (NADH, NADPH). No carcinoma renal, a reprogramação glicolítica é particularmente relevante porque se sobrepõe a alterações genéticas condutoras específicas do subtipo tumoral. A via `hsa00010` captura tanto as enzimas glicolíticas canônicas quanto enzimas da gliconeogênese e do metabolismo de aldeídos/álcoois, oferecendo uma visão integrada das alterações transcriptômicas;
5. **Conecta-se a fenótipos clinicamente relevantes**: alterações glicolíticas em tumores renais já foram associadas a prognóstico, resposta terapêutica e plasticidade metabólica;
6. **Oferece um modelo pedagógico** para alfabetização genômica: os conceitos de *upregulation*, *downregulation*, modularidade de rede e limitações do *bulk RNA-seq* podem ser ilustrados a partir de genes concretos.

### Questão de pesquisa

> Quais alterações transcriptômicas na via glicólise/gliconeogênese em KIRP podem ser identificadas a partir de dados públicos TCGA/GTEx e como esses achados podem contribuir, como exemplos formativos, para a alfabetização genômica do enfermeiro em saúde de precisão?

---

## Desenho do estudo

| Característica | Descrição |
|----------------|-----------|
| **Tipo** | Estudo computacional *in silico* com análise secundária de dados transcriptômicos públicos |
| **Fonte** | TCGA/GTEx integrados via UCSC Xena |
| **Via analisada** | KEGG `hsa00010` (glicólise/gliconeogênese) — 65 genes |
| **Grupos** | ~290 KIRP (*TCGA Kidney Papillary Cell Carcinoma*) vs. ~30 Normal (*GTEX Kidney*) |
| **Expressão diferencial** | *limma* com contraste KIRP − Normal, *eBayes*, Benjamini-Hochberg |
| **Critério DEG** | FDR < 0,05 e \|logFC\| > 1 |
| **Rede PPI** | STRING v11.5, *Homo sapiens*, *combined score* ≥ 700 |
| **Visualização** | *Volcano plot* + rede PPI (PNG, 300 dpi) |

---

## Fluxo analítico

```
┌─────────────────────────────────────────────────────────┐
│                    PIPELINE COMPLETO                      │
├───────────────┬─────────────────────────────────────────┤
│  01_download  │  Consulta KEGG hsa00010                  │
│      ↓        │  Extrai 65 símbolos gênicos humanos     │
│  02_prepare   │  Importa kidney.tsv (UCSC Xena)          │
│      ↓        │  Classifica Normal vs. KIRP              │
│               │  Converte para matriz numérica           │
│  03_diff_expr │  Ajusta modelo linear (limma)            │
│      ↓        │  Contraste KIRP - Normal                 │
│               │  Moderação eBayes + BH                   │
│               │  Classifica DEGs (FDR<0.05, |logFC|>1)   │
│  04_volcano   │  Gera volcano plot (ggplot2 + ggrepel)   │
│      ↓        │  PNG 300 dpi, genes rotulados            │
│  05_ppi       │  Mapeia DEGs → STRING IDs                │
│               │  Recupera interações (score ≥ 700)       │
│               │  Constrói rede (igraph + ggraph)         │
│               │  PNG 300 dpi, maior componente conexo    │
└───────────────┴─────────────────────────────────────────┘
```

---

## Estrutura do repositório

```text
KIRP-Glycolysis-Transcriptomics/
│
├── README.md                          ← Documentação completa
├── LICENSE                            ← MIT
├── .gitignore                         ← Bloqueia manuscript/, output/, *.Rmd, *.pdf
├── CITATION.cff                       ← Metadados de citação (GitHub)
│
├── environment/                       ← Reprodutibilidade
│   ├── sessionInfo.txt                ← Ambiente R do desenvolvimento
│   └── packages.csv                   ← Versões exatas dos pacotes
│
├── scripts/                           ← Pipeline (execução sequencial)
│   ├── 01_download_data.R             ← Download KEGG hsa00010
│   ├── 02_prepare_data.R              ← Preparo da matriz TCGA/GTEx
│   ├── 03_differential_expression.R   ← Expressão diferencial (limma)
│   ├── 04_volcano_plot.R              ← Volcano plot (ggplot2)
│   ├── 05_ppi_network.R               ← Rede PPI (STRING + igraph)
│   └── run_pipeline.R                 ← Orquestrador (roda todos em ordem)
│
├── data/                              ← Dados (não versionados no GitHub)
│   ├── raw/                           ← kidney.tsv (matriz UCSC Xena)
│   ├── processed/                     ← metadata.rds, expression_matrix.rds
│   └── metadata/                      ← Hsa_genes.csv (genes da via)
│
├── results/                           ← Resultados versionados
│   ├── differential_expression/       ← DEG_KIRP_vs_Normal.csv
│   ├── figures/                       ← Volcano.png, PPI_network.png
│   ├── ppi/                           ← STRING_mapping.csv, edges
│   └── tables/                        ← deg_summary.csv, ppi_network_summary.csv
│
└── output/                            ← Saídas auxiliares (não versionadas)
```

---

## Dependências

### R e versões dos pacotes

O *pipeline* foi desenvolvido e testado em **R 4.6.0** (Windows 11, ucrt). As versões exatas dos pacotes utilizados são:

| Pacote | Versão | Função no pipeline |
|--------|--------|--------------------|
| `KEGGREST` | 1.52.0 | Consulta à API do KEGG para obter a via `hsa00010` |
| `dplyr` | 1.2.1 | Manipulação de *data frames* e *tibbles* |
| `stringr` | 1.6.0 | Processamento de strings (extração de símbolos gênicos) |
| `tibble` | 3.3.1 | Estruturas de dados tabulares |
| `rio` | 1.3.0 | Importação/exportação de CSV, RDS e outros formatos |
| `limma` | 3.68.0 | Modelagem linear, contraste, *eBayes*, Benjamini-Hochberg |
| `ggplot2` | 4.0.3 | Gramática de gráficos (*volcano plot*) |
| `ggrepel` | 0.9.8 | Rótulos não sobrepostos no *volcano plot* |
| `STRINGdb` | 2.24.0 | Interface R para o banco STRING (PPI) |
| `igraph` | 2.3.0 | Manipulação e análise de grafos |
| `ggraph` | 2.2.2 | Visualização de grafos com gramática *ggplot2* |

### Instalação rápida

```r
install.packages(c(
  "KEGGREST", "dplyr", "stringr", "tibble", "rio",
  "ggplot2", "ggrepel", "igraph", "ggraph"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("limma", "STRINGdb"))
```

---

## Fonte dos dados

### Matriz de expressão

O arquivo `data/raw/kidney.tsv` é uma matriz de expressão gênica integrada TCGA/GTEx exportada pela plataforma **UCSC Xena** (*https://xenabrowser.net/*).

**Especificações esperadas:**

| Campo | Descrição |
|-------|-----------|
| `sample` | Identificador da amostra |
| `primary_site` | Sítio anatômico primário (*Kidney*) |
| `sample_type` | Tipo de amostra |
| `study` | Estudo de origem (TCGA ou GTEx) |
| `TCGA_GTEX_main_category` | Categoria biológica: `GTEX Kidney` ou `TCGA Kidney Papillary Cell Carcinoma` |
| Demais colunas | Expressão gênica em escala *log2(norm_count + 1)*, uma coluna por gene |

### Como obter os dados

1. Acesse *https://xenabrowser.net/datapages/*
2. Selecione o *dataset* **GDC TCGA Kidney Papillary Cell Carcinoma (KIRP)** combinado com **GTEx Kidney**
3. Faça o *download* da matriz de expressão gênica (*gene expression RNAseq — HTSeq — Counts*)
4. Renomeie o arquivo para `kidney.tsv` e coloque-o em `data/raw/`

> **Nota:** O arquivo bruto não é versionado neste repositório por questões de tamanho (~300 MB). O *pipeline* valida automaticamente a presença das colunas e categorias esperadas no script `02_prepare_data.R`.

### Dados da via metabólica

A via `hsa00010` é obtida **automaticamente** pelo script `01_download_data.R` via API do KEGG (`KEGGREST::keggGet()`). Não requer *download* manual.

---

## Como executar o pipeline

### Pré-requisitos

1. **R ≥ 4.3.0** instalado
2. Pacotes listados em [Dependências](#dependências) instalados
3. Arquivo `kidney.tsv` colocado em `data/raw/`
4. Conexão com internet (para consulta ao KEGG e STRING)

### Execução completa (recomendado)

```r
# No diretório raiz do repositório:
source("scripts/run_pipeline.R")
```

O orquestrador `run_pipeline.R`:
1. Verifica se todos os pacotes estão instalados
2. Executa os scripts 01 a 05 em sequência
3. Interrompe com mensagem de erro se algum pré-requisito falhar
4. Exibe o progresso de cada etapa no console

### Execução seletiva

Cada script pode ser executado independentemente, desde que as dependências de arquivos anteriores estejam satisfeitas:

```r
# Apenas baixar genes do KEGG:
source("scripts/01_download_data.R")

# Apenas preparar a matriz (requer kidney.tsv):
source("scripts/02_prepare_data.R")

# Apenas expressão diferencial (requer 02):
source("scripts/03_differential_expression.R")

# Apenas volcano plot (requer 03):
source("scripts/04_volcano_plot.R")

# Apenas rede PPI (requer 03 + internet):
source("scripts/05_ppi_network.R")
```

### Verificação de integridade

Após a execução, verifique se todos os arquivos de saída foram gerados:

```r
expected_outputs <- c(
  "data/metadata/Hsa_genes.csv",
  "data/processed/metadata.rds",
  "data/processed/expression_matrix.rds",
  "results/differential_expression/DEG_KIRP_vs_Normal.csv",
  "results/figures/Volcano.png",
  "results/figures/PPI_network.png",
  "results/ppi/STRING_mapping.csv",
  "results/ppi/STRING_edges_high_confidence.csv",
  "results/tables/deg_summary.csv",
  "results/tables/ppi_network_summary.csv"
)
missing <- expected_outputs[!file.exists(expected_outputs)]
if (length(missing) > 0) {
  warning("Missing outputs: ", paste(missing, collapse = ", "))
} else {
  message("All outputs generated successfully.")
}
```

---

## Descrição detalhada dos scripts

### `01_download_data.R` — Obtenção da via KEGG

Este script é responsável pela única etapa de *download* automatizado do *pipeline*.

**O que faz:**
1. Consulta a API do KEGG com `KEGGREST::keggGet("hsa00010")`
2. Extrai o campo `GENE` da resposta
3. Processa os símbolos gênicos:
   - Remove anotações entre colchetes (`[...]`)
   - Extrai apenas o primeiro símbolo antes de `;`
   - Remove duplicatas
   - Ordena alfabeticamente
4. Salva `data/metadata/Hsa_genes.csv` com 65 genes únicos

**Validações:**
- Confirma que a resposta KEGG contém o campo `GENE`
- Interrompe se a via não for encontrada

**Saída:** `data/metadata/Hsa_genes.csv`

---

### `02_prepare_data.R` — Preparo da matriz TCGA/GTEx

Este script transforma a matriz bruta da UCSC Xena em objetos prontos para análise.

**O que faz:**
1. Importa `data/raw/kidney.tsv`
2. Verifica a presença das 5 colunas de metadados obrigatórias
3. Separa metadados da matriz de expressão
4. Converte a matriz para formato numérico
5. Classifica amostras como `Normal` (GTEx) ou `KIRP` (TCGA)
6. Remove amostras que não pertencem a nenhum dos dois grupos
7. Emite *warning* se a matriz parecer estar em escala linear (inteiros), sugerindo verificação da escala *log2*
8. Salva `metadata.rds` e `expression_matrix.rds` em `data/processed/`

**Validações:**
- Arquivo `kidney.tsv` existe
- Colunas obrigatórias presentes
- Categorias `GTEX Kidney` e `TCGA Kidney Papillary Cell Carcinoma` presentes
- *Warning* de escala linear vs. logarítmica

**Saídas:**
- `data/processed/metadata.rds`
- `data/processed/expression_matrix.rds`
- `data/processed/qc_tables.rds`

---

### `03_differential_expression.R` — Expressão diferencial com limma

Este é o núcleo estatístico do *pipeline*.

**O que faz:**
1. Carrega `metadata.rds` e `expression_matrix.rds`
2. Transpõe a matriz para o formato genes × amostras
3. Constrói matriz de desenho sem intercepto (`model.matrix(~ 0 + condition)`)
4. Ajusta modelo linear com `lmFit()`
5. Define contraste `KIRP − Normal` com `makeContrasts()`
6. Aplica contraste com `contrasts.fit()`
7. Moderação empírico-bayesiana com `eBayes()` — estabiliza variâncias, especialmente importante dado o desbalanço amostral (~290 KIRP vs. ~30 Normal)
8. Extrai resultados com `topTable()` (todos os genes, ordenados por *p-value*)
9. Aplica correção de Benjamini-Hochberg (FDR)
10. Classifica genes: `Up` (FDR < 0,05 e logFC > 1), `Down` (FDR < 0,05 e logFC < −1), `NS` (não significativo)

**Por que matriz sem intercepto?**
A parametrização `~ 0 + condition` estima a média de expressão de cada grupo independentemente. Isso facilita a interpretação biológica dos coeficientes (média do grupo Normal, média do grupo KIRP) antes da aplicação do contraste. A alternativa com intercepto produziria coeficientes como "intercepto + efeito KIRP", que são menos intuitivos para este desenho simples de dois grupos. Ambas as parametrizações produzem resultados de contraste idênticos.

**Por que |logFC| > 1?**
O limiar de |logFC| > 1 equivale a um *fold change* mínimo de 2× (2¹ = 2). Este limiar é convenção amplamente adotada em estudos transcriptômicos de vias metabólicas (ver, por exemplo, Lv et al., 2021; Bao et al., 2020), equilibrando relevância estatística (FDR < 0,05) com magnitude biológica do efeito. O *pipeline* preserva os valores de logFC e FDR para todos os genes na tabela completa, permitindo que o usuário aplique limiares alternativos sem reexecutar a modelagem.

**Saídas:**
- `results/differential_expression/DEG_KIRP_vs_Normal.csv` (todos os genes com logFC, FDR, classificação)
- `results/tables/deg_summary.csv` (contagem: Up, Down, NS)

---

### `04_volcano_plot.R` — Visualização da expressão diferencial

**O que faz:**
1. Carrega `DEG_KIRP_vs_Normal.csv`
2. Aplica os mesmos limiares (FDR < 0,05, |logFC| > 1) para classificação visual
3. Constrói *volcano plot* com `ggplot2`:
   - Eixo x: logFC
   - Eixo y: `-log10(FDR)`
   - Linhas tracejadas nos limiares de decisão
   - Cores: azul (*up*), roxo (*down*), cinza (NS)
4. Rotula todos os DEGs com `ggrepel` (evita sobreposição)
5. Exporta PNG com 300 dpi

**Parâmetros ajustáveis:** `lfc_cutoff` e `fdr_cutoff` no início do script.

**Saída:** `results/figures/Volcano.png` (7×6 polegadas, 300 dpi)

---

### `05_ppi_network.R` — Rede de interação proteína-proteína

**O que faz:**
1. Carrega `DEG_KIRP_vs_Normal.csv`
2. Separa genes *upregulated* e *downregulated*
3. Conecta-se ao STRING (v11.5, *Homo sapiens*, taxID 9606)
4. Mapeia símbolos gênicos → identificadores STRING
5. Recupera todas as interações entre os DEGs mapeados
6. Filtra arestas com `combined_score ≥ 700` (alta confiança)
7. Constrói grafo não direcionado com `igraph`
8. Extrai o maior componente conexo
9. Visualiza com `ggraph` (layout Fruchterman-Reingold, `set.seed(1)` para reprodutibilidade):
   - Espessura das arestas proporcional ao *score*
   - Cores dos nós por direção de regulação
   - Rótulos com `ggrepel`
10. Exporta PNG com 300 dpi
11. Salva tabelas de mapeamento, arestas e sumário da rede

**Parâmetros ajustáveis:**
- `score_cutoff` (atual: 700)
- `version` do STRING (atual: 11.5)
- `species` (atual: 9606)

**Saídas:**
- `results/figures/PPI_network.png` (9×7 polegadas, 300 dpi)
- `results/ppi/STRING_mapping.csv`
- `results/ppi/STRING_edges_high_confidence.csv`
- `results/tables/ppi_network_summary.csv`

---

## Saídas geradas

### Figuras

| Arquivo | Descrição | Resolução |
|---------|-----------|-----------|
| `results/figures/Volcano.png` | Genes *up/down*/NS com limiares FDR e logFC | 300 dpi |
| `results/figures/PPI_network.png` | Maior componente conexo, colorido por regulação | 300 dpi |

### Tabelas

| Arquivo | Conteúdo |
|---------|----------|
| `data/metadata/Hsa_genes.csv` | 65 genes da via `hsa00010` |
| `results/differential_expression/DEG_KIRP_vs_Normal.csv` | Todos os genes com logFC, AveExpr, t, P.Value, adj.P.Val, B, regulation |
| `results/tables/deg_summary.csv` | Contagem: Up × Down × NS |
| `results/ppi/STRING_mapping.csv` | Mapeamento gene_symbol → STRING_id |
| `results/ppi/STRING_edges_high_confidence.csv` | Arestas com combined_score ≥ 700 |
| `results/tables/ppi_network_summary.csv` | Métricas: nós, arestas, componente conexo |

---

## Escopo e limitações analíticas

### O que este pipeline **faz**

- Análise de expressão diferencial com *limma* (FDR + logFC)
- *Volcano plot* com genes rotulados
- Rede PPI com STRING (alta confiança)
- Documentação completa de cada etapa

### O que este pipeline **NÃO faz**

- Análise de sobrevida (Kaplan-Meier, Cox)
- Enriquecimento funcional (GO, KEGG *pathway* enrichment)
- Análise de *splicing* alternativo ou isoformas
- PPI tridimensional ou análise estrutural de proteínas
- UpSet plot, *dotplot* ou *heatmap*
- Validação experimental (proteica, funcional)
- Normalização ou controle de qualidade *de novo* dos dados brutos (a matriz UCSC Xena já é pré-processada)

### Limitações metodológicas

1. **Nível transcriptômico apenas:** RNA não equivale a proteína funcional. Modificações pós-transcricionais, *splicing* alternativo, regulação traducional e atividade enzimática não são capturadas.
2. ***Bulk RNA-seq*:** a expressão medida reflete a média de populações celulares heterogêneas (células tumorais, estroma, imunes, vasculares). Variações célula a célula não são detectadas.
3. **Desbalanço amostral:** ~290 KIRP vs. ~30 Normal. O *limma* mitiga esse problema com moderação empírico-bayesiana, mas a estimativa de variância no grupo normal é menos precisa.
4. **Rede PPI *in silico*:** interações do STRING são preditas ou curadas a partir de evidências heterogêneas. A presença de uma aresta não demonstra interação física *in vivo* no contexto do KIRP.
5. **Generalização entre subtipos:** a maior parte da literatura de reprogramação metabólica em câncer renal deriva de KIRC (células claras). Os mecanismos moleculares podem diferir no KIRP.

---

## Reprodutibilidade

### O que garante a reprodutibilidade

- **Scripts independentes e sequenciais:** cada etapa é autocontida, com validações de entrada
- **Versões de pacotes documentadas:** `environment/packages.csv` e `sessionInfo.txt`
- **Semente fixa:** `set.seed(1)` no layout da rede PPI
- **Parâmetros explícitos:** todos os limiares e opções declarados no código
- **Validações de integridade:** cada script verifica se os arquivos de entrada existem antes de prosseguir

### Limitações à reprodutibilidade

- O arquivo `kidney.tsv` (~300 MB) não é versionado — requer *download* manual da UCSC Xena
- A consulta ao STRING depende de conexão com internet e da disponibilidade do servidor
- Pequenas variações numéricas podem ocorrer entre versões de pacotes (embora *limma*, *ggplot2* e *igraph* sejam estáveis)

---

## Declaração de uso de Inteligência Artificial

Em conformidade com a **Portaria CNPq nº 2.664, de 22 de dezembro de 2026**, declara-se que ferramentas de inteligência artificial generativa foram utilizadas como recursos de apoio durante a elaboração deste projeto. Especificamente, foram empregados:

- **ChatGPT-5.5** (OpenAI)
- **Codex** (OpenAI)
- **DeepSeek-v4-pro** (Hangzhou DeepSeek Artificial Intelligence)

nas seguintes etapas: (i) revisão gramatical, estilística e de coesão textual da documentação; (ii) sugestão de estrutura e organização do repositório; (iii) verificação de consistência entre referências e conteúdo; (iv) formatação de referências conforme o estilo Vancouver; e (v) elaboração de documentação técnica.

Em nenhum momento as ferramentas de IA substituíram o julgamento científico dos autores, a análise crítica dos resultados, a interpretação dos dados transcriptômicos ou a responsabilidade pelo conteúdo final. Todos os autores revisaram, editaram e aprovam integralmente o conteúdo, assumindo total responsabilidade intelectual pelo repositório.

---

## Licença

Este projeto está licenciado sob a **MIT License**. Consulte o arquivo [LICENSE](LICENSE) para detalhes.

---

## Contato

**Ryan de Paulo Santos**
ryan.paulo@gsuite.iff.edu.br

Instituto Federal de Educação, Ciência e Tecnologia Fluminense — *Campus* Campos Guarus

---

## Como citar

Use o arquivo [CITATION.cff](CITATION.cff) para citação automática pelo GitHub.

**Formato sugerido (Vancouver):**

> Loureiro KC, Santos RP, Freitas LMD, Ivine, Pecly ME. KIRP Glycolysis Transcriptomics [software]. GitHub; 2026 [citado 2026 jul 8]. Disponível em: https://github.com/santosry/KIRP-Glycolysis-Transcriptomics.
