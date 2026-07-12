# Central Carbon Metabolism Transcriptomics in KIRP

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R 4.6.0](https://img.shields.io/badge/R-4.6.0-276DC3)](https://www.r-project.org/)
[![v3.1.0](https://img.shields.io/badge/version-3.1.0-green)]()
[![Python 3](https://img.shields.io/badge/Python-3.x-3776AB)](https://www.python.org/)
[![CI](https://img.shields.io/badge/CI-validated-brightgreen)]()

**Perfil transcriptômico dos genes do metabolismo central do carbono no carcinoma renal papilar (KIRP): análise pareada tumor-tecido adjacente, sensibilidade à escolha do tecido de referência e implicações metodológicas.**

*Transcriptomic profile of central carbon metabolism genes in papillary renal cell carcinoma (KIRP): paired tumor-adjacent analysis, sensitivity to reference tissue choice, and methodological implications.*

---

## Índice

1. [Questão Científica](#questão-científica)
2. [Vias Analisadas](#vias-analisadas)
3. [Principais Achados](#principais-achados)
4. [Desenho dos Comparadores](#desenho-dos-comparadores)
5. [Dados](#dados)
6. [Estrutura do Repositório](#estrutura-do-repositório)
7. [Reprodução](#reprodução)
8. [Instruções para Agentes de IA](#instruções-para-agentes-de-ia)
9. [Visualizações 3D Interativas](#visualizações-3d-interativas)
10. [Limitações](#limitações)
11. [Declaração de Uso de Inteligência Artificial](#declaração-de-uso-de-inteligência-artificial)
12. [Como Citar](#como-citar)
13. [Referências](#referências)
14. [Licença](#licença)

---

## Questão Científica

Qual é o perfil transcriptômico dos genes do metabolismo central do carbono (glicólise/gliconeogênese, via das pentoses fosfato e ciclo do ácido cítrico) no carcinoma renal papilar (KIRP) quando analisado com tecido normal adjacente pareado como comparador primário, e qual a sensibilidade dos resultados à escolha do tecido de referência?

*What is the transcriptomic profile of central carbon metabolism genes (glycolysis/gluconeogenesis, pentose phosphate pathway, citrate cycle) in KIRP when analyzed with paired tumor-adjacent normal tissue as the primary comparator, and how sensitive are the results to the choice of reference tissue?*

---

## Vias Analisadas

| KEGG ID | Via / Pathway | Genes no KEGG | Genes na Matriz |
|---------|---------------|:-------------:|:---------------:|
| hsa00010 | Glicólise / Gliconeogênese (Glycolysis / Gluconeogenesis) | 67 | 64 |
| hsa00030 | Via das Pentoses Fosfato (Pentose Phosphate Pathway) | 31 | 30 |
| hsa00020 | Ciclo do Ácido Cítrico — TCA (Citrate Cycle) | 30 | 29 |
| **União (únicos)** | **Metabolismo Central do Carbono** | **110** | **106** |

**Genes ausentes da matriz (n = 4):** *G6PC1*, *PRPS1L1*, *RPEL1*, *SUCLA2* — removidos durante a filtragem de baixa expressão (expressão ≤ 1 em > 90% das amostras). Verificação de sensibilidade: a remoção destes genes não afetou quaisquer conclusões.

## Por que estas três vias?

### Glicólise / Gliconeogênese (hsa00010) — o eixo central

A via glicolítica é o ponto de partida do metabolismo central do carbono: converte glicose em piruvato, gerando ATP e intermediários biossintéticos que alimentam todas as demais vias metabólicas. A gliconeogênese, via reversa, é crítica no rim pois o túbulo proximal renal é um dos principais sítios gliconeogênicos do organismo humano. No câncer, o efeito Warburg — aumento da captação de glicose e fermentação a lactato mesmo na presença de oxigênio — é uma das alterações metabólicas mais bem documentadas (VANDER HEIDEN et al., 2009). Estudar esta via no KIRP é essencial por duas razões: (i) o rim normal expressa enzimas gliconeogênicas (*FBP1*, *PCK1*, *ALDOB*) em altíssimos níveis no túbulo proximal, e sua perda no tumor serve como marcador de desdiferenciação epitelial; (ii) isoenzimas glicolíticas como *HK2* e *PKM2* são alvos terapêuticos em múltiplos tipos de câncer e sua expressão no KIRP não havia sido sistematicamente caracterizada em análise pareada. Com 64 genes na matriz, é a maior das três vias e fornece o contexto transcriptômico principal.

### Via das Pentoses Fosfato (hsa00030) — defesa redox e biossíntese

A PPP ramifica-se da glicólise na glicose-6-fosfato e desempenha duas funções essenciais para a proliferação celular: (i) produção de NADPH, o principal redutor celular que sustenta a defesa antioxidante (glutationa, tiorredoxina) e a biossíntese de lipídeos; (ii) produção de ribose-5-fosfato, precursor da síntese de nucleotídeos (DNA, RNA). Em carcinomas renais, a PPP é particularmente relevante porque a enzima limitante G6PD é regulada pela via NRF2 — fator de transcrição frequentemente ativado em KIRP com deficiência de FH (KOPPULA et al., 2020; PATRA; HAY, 2014). Além disso, a produção de NADPH é crítica para a homeostase redox em um órgão exposto a alto estresse oxidativo como o rim. A inclusão desta via permite testar se há evidência transcriptômica de ativação coordenada do braço oxidativo da PPP no KIRP. Com 30 genes na matriz, a PPP complementa a via glicolítica com enzimas que determinam o destino do carbono entre oxidação (geração de NADPH) e biossíntese (geração de ribose).

### Ciclo do Ácido Cítrico — TCA (hsa00020) — encruzilhada metabólica e tumorigênese

O ciclo do ácido cítrico (ciclo de Krebs) é a encruzilhada central do metabolismo oxidativo: oxida acetil-CoA a CO₂, gerando equivalentes redutores (NADH, FADH₂) para a cadeia respiratória, e fornece intermediários para biossíntese de aminoácidos, lipídeos e heme. No contexto específico do KIRP, o TCA tem relevância única: mutações no gene *FH* (fumarato hidratase), que codifica uma enzima do ciclo, são a base molecular do subtipo hereditário de KIRP com deficiência de FH (HLRCC). O acúmulo de fumarato promove succinação de PTEN, conectando diretamente o TCA à tumorigênese neste subtipo (GE et al., 2022). Além disso, vários genes do ciclo (CS, ACO2, IDH2/3, OGDH, SDHA-D, MDH2) codificam enzimas mitocondriais cuja expressão pode refletir o conteúdo mitocondrial e o estado de diferenciação celular. A inclusão sistemática de todos os 29 genes do TCA disponíveis na matriz permite: (i) verificar se há alteração coordenada do ciclo como conjunto gênico — hipótese testada formalmente via *camera* — e (ii) contextualizar os achados de *FH*, *SDHs* e *IDHs* cujas mutações são recorrentes em neoplasias renais (SANCHEZ; SIMON, 2024).

### Integração das três vias

As três vias não operam isoladamente: a glicólise alimenta a PPP (via glicose-6-fosfato) e o TCA (via piruvato → acetil-CoA); a PPP fornece NADPH para biossíntese redutora e ribose-5-fosfato para nucleotídeos; o TCA fornece citrato para lipogênese e intermediários para aminoácidos. Analisar as três vias em conjunto — 106 genes extraídos *a posteriori* de um modelo transcriptômico único sobre 31.633 genes — permite capturar o estado coordenado do metabolismo central do carbono sem o viés de seleção que ocorreria ao restringir previamente o universo gênico apenas a essas vias (circularidade metodológica). A análise conjunta também revela genes compartilhados entre vias (7 dos 35 DEGs pertencem a duas ou três vias), evidenciando pontos de conexão metabólica potencialmente relevantes.

---

## Principais Achados (v3.1.0)

### Análise pareada primária (32 pares tumor-adjacente KIRP)

- **35/106 genes (33,0%)** diferencialmente expressos (|log₂FC| > 1, FDR < 0,05): 9 aumentados, 26 diminuídos
- 7 genes compartilhados entre duas ou três vias metabólicas
- **Glicólise (64 genes):** 26 DEGs (7 aumentados, 19 diminuídos)
- **Via das Pentoses Fosfato (30 genes):** 12 DEGs (3 aumentados, 9 diminuídos)
- **TCA (29 genes):** 4 DEGs (todos diminuídos)

### Genes de maior magnitude

| Gene | log₂FC | IC 95% | FDR | Regulação |
|------|:------:|--------|:---:|:---------:|
| *ALDOB* (aldolase B) | -8,66 | [-9,96; -7,36] | 4,3×10⁻¹³ | ↓ Diminuído |
| *ADH1C* | -6,01 | [-7,06; -4,95] | 1,9×10⁻¹⁵ | ↓ Diminuído |
| *PCK1* (PEPCK) | -5,11 | [-6,53; -3,69] | 1,7×10⁻⁹ | ↓ Diminuído |
| *HK2* (hexoquinase 2) | +3,34 | [2,51; 4,17] | 1,9×10⁻¹² | ↑ Aumentado |
| *FBP1* (FBPase 1) | -3,07 | [-3,90; -2,23] | 1,1×10⁻⁸ | ↓ Diminuído |
| *G6PD* (G6PDH) | +1,49 | [0,89; 2,09] | 1,8×10⁻⁶ | ↑ Aumentado |
| *TKT* (transcetolase) | +1,60 | [1,19; 2,00] | 8,1×10⁻¹² | ↑ Aumentado |

### Teste de conjuntos gênicos (camera)

- **hsa00020 (TCA):** deslocamento coordenado no sentido de redução — FDR = 0,0012 (*inter.gene.cor = 0,01*)
- **hsa00010 (Glicólise):** FDR = 0,078 (sem evidência estatística após correção)
- **hsa00030 (PPP):** FDR = 0,800 (sem deslocamento detectável)
- **Sensibilidade:** resultados não se mantêm com correlação estimada por conjunto (*inter.gene.cor = NA*)

### Concordância entre comparadores

| Métrica | Pareado vs. TCGA-KIRP expandida | Pareado vs. GTEx |
|---------|:-------------------------------:|:----------------:|
| CCC de Lin | 0,974 | 0,795 |
| MAE (log₂FC) | 0,27 | 0,94 |
| Viés (d = pareado − comparador) | 0,00 | −0,79 |
| Limites de concordância 95% | [−0,53; 0,54] | [−2,24; 0,67] |
| Concordância direcional | 98,7% (74/75) | 60,2% (59/98) |
| Genes discordantes (todos) | 4 | 45 |

> ⚠️ **Interpretação:** GTEx como referência externa introduziu viés sistemático de −0,79, superestimando os log₂FC e revertendo a direção em 45 dos 106 genes. As duas análises com tecido adjacente KIRP compartilham os mesmos 32 controles; a alta concordância **não** representa validação independente.

---

## Desenho dos Comparadores

| Comparador | Delineamento | n | Papel |
|------------|-------------|:-:|-------|
| **Pareado KIRP** | 32 pares tumor-adjacente, modelo: ~ patient + condition | 64 | **Primário** |
| KIRP vs. Adjacente Normal | 288 tumores vs. 32 KIRP-adjacentes (não pareado) | 320 | Secundário (mesma coorte) |
| KIRP vs. Todos TCGA Normal | 288 tumores vs. 129 normais TCGA | 417 | Exploratório (projetos mistos) |
| KIRP vs. GTEx Normal | Cross-coorte, estruturalmente confundido | 316 | **Sensibilidade apenas** |

---

## Dados

### Fonte

| Propriedade | Valor |
|-------------|-------|
| **Plataforma** | [UCSC Xena](https://xenabrowser.net/) |
| **Dataset** | `TcgaTargetGtex_RSEM_Hugo_norm_count` |
| **Hub** | `toil-xena-hub` (pipeline Toil) |
| **Amostras** | 445 (288 KIRP, 28 GTEx Normal, 32 KIRP-adjacente, 97 outros TCGA) |
| **Genes (brutos)** | 58.581 |
| **Genes (após filtragem)** | 31.633 (26.948 removidos por baixa expressão) |
| **Escala** | log₂(norm_count + 1) — pré-transformada pelo UCSC Xena |
| **Faixa de expressão** | 0 a 20,1 (mediana = 11,15) |
| **Acesso** | 11 de julho de 2026 |
| **SHA256 (kidney.tsv)** | `29154696504ce365e681d9c319fe352a6c84c5ae87798c8f1cce3c11159f7ea2` |
| **SHA256 (kidney_transcriptome.tsv)** | consulte `results/v3/checksums_sha256.txt` |
| **Tamanho (transcriptoma)** | ~178 MB (TSV); ~70 MB (gzip) |

### Como obter os dados

#### Opção 1: Baixar o arquivo kidney.tsv (genes das vias metabólicas, ~220 KB)

Este arquivo está versionado neste repositório em `data/raw/kidney.tsv`. Para baixá-lo diretamente do UCSC Xena:

1. Acesse: [https://xenabrowser.net/datapages/](https://xenabrowser.net/datapages/)
2. Selecione o dataset: **GDC TCGA Kidney Papillary Cell Carcinoma (KIRP)** combinado com **GTEx Kidney**
3. Escolha a matriz: **gene expression RNAseq — RSEM norm_count (log₂ transformed)**
4. Clique em "Download" e salve como `data/raw/kidney.tsv`

Alternativamente, utilize o script Python incluso:

```bash
# Baixa o dataset kidney (painel gênico reduzido)
# Este arquivo é usado pelo download_full_transcriptome.py como referência de amostras
python3 scripts/download_full_matrix.py
```

#### Opção 2: Baixar o transcriptoma completo (~178 MB, ~70 MB zipado)

```bash
# Requer Python 3 e ~2 GB de espaço temporário durante o download
# Baixa todos os 58.581 genes, filtra para amostras renais
python3 scripts/download_full_transcriptome.py
```

**Verificação de integridade:**

```bash
sha256sum data/raw/kidney.tsv
# Esperado: 29154696504ce365e681d9c319fe352a6c84c5ae87798c8f1cce3c11159f7ea2

sha256sum data/raw/kidney_transcriptome.tsv
# Consulte results/v3/checksums_sha256.txt para o hash esperado
```

---

## Estrutura do Repositório

```
KIRP-Glycolysis-Transcriptomics/
│
├── data/
│   ├── raw/                            # Dados brutos (alguns versionados)
│   │   ├── kidney.tsv                  # Painel gênico reduzido (~220 KB, versionado)
│   │   └── kidney_transcriptome.tsv    # Transcriptoma completo (~178 MB, NÃO versionado)
│   ├── processed/                      # Dados processados (RDS)
│   │   ├── expression_hsa00010.rds     # Matriz glicólise
│   │   ├── expression_hsa00020.rds     # Matriz TCA
│   │   ├── expression_hsa00030.rds     # Matriz PPP
│   │   └── expression_matrix_full.rds  # Matriz completa (amostras x genes)
│   ├── metadata/                       # Metadados de genes
│   └── provenance/                     # Registro de proveniência
│
├── scripts/
│   ├── pipeline_v3.R                   # Pipeline principal v3 (análise pareada)
│   ├── pipeline_v3_addendum.R          # Adendo: concordância, modelo não-pareado-64
│   ├── 16_3d_visualizations.R          # Volcano plots 3D interativos (HTML)
│   ├── 16b_ppi_3d_correlation.R        # Rede de coexpressão 3D (HTML)
│   ├── download_full_transcriptome.py  # Download do transcriptoma completo
│   ├── download_full_matrix.py         # Download do painel reduzido
│   ├── 00_environment.R                # Verificação de dependências
│   └── *.R                             # Scripts auxiliares (v2, auditoria, etc.)
│
├── results/
│   ├── v3/                             # Resultados da versão 3
│   │   ├── tables/                     # DEGs, tabelas suplementares, concordância
│   │   │   ├── Supplementary_Table_S1.csv  # Todos os 106 genes
│   │   │   ├── DEG_hsa00010.csv        # DEGs da via glicolítica
│   │   │   ├── DEG_hsa00020.csv        # DEGs do TCA
│   │   │   ├── DEG_hsa00030.csv        # DEGs da PPP
│   │   │   ├── discordant_genes.csv    # Genes com direção discordante (Paired vs GTEx)
│   │   │   ├── camera_gene_sets.csv    # Resultados do teste camera
│   │   │   ├── ppi_3d_centrality.csv   # Centralidade da rede 3D
│   │   │   └── ppi_3d_summary.csv      # Resumo da rede 3D
│   │   ├── figures/
│   │   │   ├── Volcano_hsa00010.png / _3D.html   # Volcano glicólise (PNG + 3D HTML)
│   │   │   ├── Volcano_hsa00020.png / _3D.html   # Volcano TCA
│   │   │   ├── Volcano_hsa00030.png / _3D.html   # Volcano PPP
│   │   │   ├── PPI_network_3D.html              # Rede PPI/coexpressão 3D
│   │   │   ├── BlandAltman_Paired_vs_GTEx.png   # Bland-Altman
│   │   │   ├── PCA_transcriptome.png            # PCA transcriptoma completo
│   │   │   ├── plotSA_paired.png                # Diagnóstico limma
│   │   │   └── Paired_*.png                     # Gráficos pareados por gene
│   │   ├── supplementary/              # Material suplementar (S1-S8)
│   │   ├── sessionInfo.txt             # Informação da sessão R
│   │   └── checksums_sha256.txt        # Checksums de todos os outputs
│   ├── differential_expression/        # Resultados v2 (legado)
│   ├── enrichment/                     # Enriquecimento funcional
│   ├── figures/                        # Figuras v2 (legado)
│   └── tables/                         # Tabelas v2 (legado)
│
├── tests/
│   └── testthat/
│       └── test_pipeline.R             # Testes automatizados (30 asserções)
│
├── environment/
│   ├── packages.csv                    # Versões dos pacotes R
│   └── sessionInfo.txt                 # Sessão R completa
│
├── renv.lock                           # Lockfile do renv
├── renv/                               # Ambiente renv (bibliotecas versionadas)
├── .github/workflows/ci.yml            # CI (GitHub Actions)
├── .gitignore
├── .Rprofile
├── LICENSE                             # MIT
├── CITATION.cff                        # Metadados de citação
├── VERSION                             # v3.1.0
├── README.md                           # Este arquivo
├── REPRODUCIBILITY.md                  # Guia de reprodutibilidade
├── DATA_PROVENANCE.md                  # Proveniência dos dados
└── manuscrito_kirp_*.Rmd / *.pdf       # Manuscrito (NÃO versionado)
```

---

## Reprodução

### ⚠️ Pré-requisitos

- **R** ≥ 4.6.0 com Bioconductor 3.23
- **Python** ≥ 3.8 (apenas para download dos dados)
- **Pacotes R:** `limma`, `dplyr`, `rio`, `ggplot2`, `ggrepel`, `pheatmap`, `igraph`, `plotly`, `htmlwidgets`, `clusterProfiler`, `org.Hs.eg.db` e dependências (lista completa em `environment/packages.csv`)
- **Pacotes Python:** `urllib` (biblioteca padrão; nenhuma instalação adicional necessária)
- **Conexão com internet:** necessária para download dos dados e anotações KEGG
- **Espaço em disco:** ~2 GB (durante download e processamento); ~250 MB (resultados finais)
- **Memória RAM:** ≥ 8 GB recomendado para o transcriptoma completo

### Ambiente R (renv)

```r
# Instalar renv se necessário
install.packages("renv")

# Restaurar ambiente exato
renv::restore()
```

### Pipeline completo (versão 3)

```bash
# 1. Clonar o repositório
git clone https://github.com/santosry/KIRP-Glycolysis-Transcriptomics.git
cd KIRP-Glycolysis-Transcriptomics

# 2. Baixar o transcriptoma completo (Python)
python3 scripts/download_full_transcriptome.py

# 3. (Opcional) Verificar checksums
sha256sum -c results/v3/checksums_sha256.txt

# 4. Executar pipeline v3 (análise pareada primária)
Rscript scripts/pipeline_v3.R

# 5. Executar adendo (análises de sensibilidade, concordância)
Rscript scripts/pipeline_v3_addendum.R

# 6. Gerar visualizações 3D interativas
Rscript scripts/16_3d_visualizations.R
Rscript scripts/16b_ppi_3d_correlation.R

# 7. Executar testes automatizados
R -e "testthat::test_dir('tests/testthat')"
```

### Passo a passo mínimo (apenas replicação dos resultados principais)

```r
# Dentro do R:
source("scripts/00_environment.R")       # Verificar dependências
source("scripts/pipeline_v3.R")          # Análise pareada, volcano, concordância
```

### Verificação de integridade dos resultados

```bash
# Verificar todos os checksums
cd results/v3
sha256sum -c checksums_sha256.txt
```

### Componentes determinísticos

| Componente | Reprodutibilidade |
|------------|:-----------------:|
| limma eBayes | ✅ Determinístico |
| PCA (prcomp) | ✅ Determinístico |
| Correlação de Pearson | ✅ Determinístico |
| Teste exato de Fisher (ORA) | ✅ Determinístico |
| Teste camera | ✅ Determinístico |
| Layout Fruchterman-Reingold 3D | ✅ Determinístico (seed fixa) |
| STRING PPI | ⚠️ Requer internet; API pode mudar |
| Anotações KEGG (KEGGREST) | ⚠️ Requer internet; release pode mudar |

---

## Instruções para Agentes de IA

> **Esta seção é destinada a agentes de inteligência artificial (Codex, Claude, ChatGPT, DeepSeek, Copilot, etc.) que estejam executando este repositório. Leia atentamente antes de executar qualquer comando.**

### Regras gerais

1. **Ambiente:** Todos os scripts devem ser executados a partir do diretório raiz do repositório (`KIRP-Glycolysis-Transcriptomics/`). O diretório de trabalho do R **deve** ser este.

2. **RENV:** Este projeto utiliza `renv` para gerenciamento de dependências. Antes de executar qualquer script R, verifique se o ambiente está sincronizado:
   ```r
   renv::status()
   renv::restore()  # se necessário
   ```
   Caso `renv` não esteja disponível, instale os pacotes listados em `environment/packages.csv` manualmente. **Não modifique o `renv.lock` a menos que explicitamente solicitado.**

3. **Ordem de execução obrigatória:** Os scripts devem ser executados nesta ordem:
   ```
   1. scripts/download_full_transcriptome.py  (Python — baixa dados)
   2. scripts/00_environment.R                 (verifica dependências)
   3. scripts/pipeline_v3.R                    (análise principal)
   4. scripts/pipeline_v3_addendum.R           (análises de sensibilidade)
   5. scripts/16_3d_visualizations.R           (volcano 3D HTML)
   6. scripts/16b_ppi_3d_correlation.R         (rede PPI 3D HTML)
   ```

4. **Dados de entrada obrigatórios:**
   - `data/raw/kidney.tsv` — **já está versionado** (~220 KB). Contém os IDs das amostras renais. **Não o sobrescreva**.
   - `data/raw/kidney_transcriptome.tsv` — **NÃO está versionado** (~178 MB). Deve ser baixado pelo script `download_full_transcriptome.py`.

5. **Timeout:** O download do transcriptoma completo (`download_full_transcriptome.py`) pode levar **10-20 minutos** dependendo da conexão. Configure timeouts adequados (≥ 1200 segundos).

6. **Pandoc não é necessário** para gerar os HTMLs 3D (os widgets são salvos com `selfcontained = FALSE` e carregam plotly.js de um diretório local). Os arquivos `*_files/` acompanhantes **devem ser mantidos junto aos HTMLs**.

7. **Não execute scripts com prefixo `10_` (`10_string_network.R`, `10b_string_network_retry.R`) a menos que explicitamente solicitado** — eles dependem da API STRING, requerem ~500 MB de download e são instáveis. A rede PPI 3D atual (`16b_ppi_3d_correlation.R`) usa correlação de Pearson (offline, reprodutível).

8. **Scripts legados (v2):** Os scripts numerados de `01` a `15` são da versão 2 do pipeline (apenas glicólise, 66 genes). **Use os scripts `pipeline_v3.R` e `16_*.R` para a versão atual (106 genes, 3 vias).**

9. **Não modifique** os arquivos em `results/v3/` manualmente. Todos são gerados automaticamente pelo pipeline.

10. **Manuscrito:** O manuscrito (`.Rmd`, `.pdf`) **não deve ser versionado** nem publicado. Se encontrar arquivos de manuscrito rastreados pelo Git, remova-os com `git rm --cached`.

11. **Ao reportar resultados:** Sempre cite os números exatos como aparecem nos arquivos de saída (`results/v3/tables/`), não os do README (que podem estar desatualizados). Confira `results/v3/sessionInfo.txt` para o ambiente exato.

12. **Checksums:** Após gerar resultados, verifique a integridade com `sha256sum -c results/v3/checksums_sha256.txt`. Atualize o arquivo de checksums se novos resultados forem gerados.

13. **Testes:** Execute `testthat::test_dir("tests/testthat")` após cada alteração no pipeline. Todos os 30 testes devem passar.

### Checklist de verificação para agentes de IA

- [ ] `data/raw/kidney.tsv` existe e tem ~220 KB
- [ ] `data/raw/kidney_transcriptome.tsv` existe e tem ~178 MB (ou foi baixado com sucesso)
- [ ] `renv::status()` retorna "The project is synchronized"
- [ ] `source("scripts/00_environment.R")` não produz erros
- [ ] `pipeline_v3.R` executa até o fim sem erros
- [ ] `pipeline_v3_addendum.R` executa até o fim sem erros
- [ ] Arquivos `results/v3/tables/Supplementary_Table_S1.csv` e `DEG_hsa*.csv` existem
- [ ] Arquivos `results/v3/figures/Volcano_*_3D.html` existem
- [ ] Arquivo `results/v3/figures/PPI_network_3D.html` existe
- [ ] `testthat::test_dir("tests/testthat")` — 30/30 passam
- [ ] Nenhum arquivo de manuscrito está rastreado pelo Git

---

## Visualizações 3D Interativas

Os arquivos HTML na pasta `results/v3/figures/` contêm visualizações 3D interativas que podem ser abertas em qualquer navegador moderno:

### Volcano Plots 3D

- **`Volcano_hsa00010_3D.html`** — Volcano 3D da via glicolítica (64 genes)
- **`Volcano_hsa00020_3D.html`** — Volcano 3D do ciclo do ácido cítrico (29 genes)
- **`Volcano_hsa00030_3D.html`** — Volcano 3D da via das pentoses fosfato (30 genes)

**Eixos:**
- **X:** log₂(Fold Change) — KIRP vs. Normal Adjacente
- **Y:** −log₁₀(FDR) — significância estatística
- **Z:** AveExpr — expressão média do gene

**Cores:** Azul (▲) = aumentado; Roxo (▼) = diminuído; Cinza = não significativo

**Interações:**
- **Arrastar:** rotaciona a visualização
- **Scroll:** zoom in/out
- **Clique duplo:** reseta a visão
- **Hover:** informações detalhadas do gene (log₂FC, FDR, IC 95%, expressão média)

### Rede PPI / Coexpressão 3D

- **`PPI_network_3D.html`** — Rede de coexpressão gênica dos 109 genes do metabolismo central do carbono

**Método:** Correlação de Pearson entre todos os pares de genes (|r| > 0,6), layout Fruchterman-Reingold 3D com seed fixa.

**Características:**
- **Nós proporcionais ao degree** (número de conexões)
- **Cores:** Azul = aumentado; Roxo = diminuído; Cinza = NS
- **Arestas:** correlações |r| > 0,6 (935 arestas, 85 nós no componente gigante)

> **Nota:** Os arquivos `*_files/` que acompanham os HTMLs contêm as bibliotecas JavaScript (plotly.js) e **devem ser mantidos junto com os HTMLs**. Não os delete nem os mova separadamente.

---

## Limitações

### Limitações do desenho experimental

1. **Tamanho amostral da análise pareada:** Apenas 32 pares tumor-adjacente. Pequenos efeitos (|log₂FC| < 1) podem não ser detectados por falta de poder estatístico. Efeitos com magnitude menor que 1,5 podem ter estimativas instáveis.

2. **Bulk RNA-seq:** Os dados são de RNA-seq de tecido total (*bulk*), que não distingue tipos celulares (células tumorais, imunes, estromais, tubulares normais). A redução de genes como *ALDOB*, *FBP1* e *PCK1* pode refletir perda de células epiteliais tubulares diferenciadas por substituição com células tumorais e estromais, e não reprogramação metabólica direcionada. **Sem estimativas de pureza tumoral ou deconvolução celular, esta hipótese concorrente não pode ser excluída.**

3. **Ausência de validação externa:** Não há coorte de validação independente. Os resultados são específicos da coorte TCGA-KIRP analisada.

### Limitações da medição

4. **mRNA ≠ proteína ≠ atividade enzimática ≠ fluxo metabólico:** A abundância de transcritos não equivale a proteína funcional, atividade enzimática ou fluxo através das vias metabólicas. Inferências sobre atividade metabólica a partir de dados transcriptômicos são indiretas e requerem validação ortogonal (proteômica, metabolômica).

5. **Dados pré-transformados:** A matriz do UCSC Xena é fornecida em escala log₂(norm_count + 1). A transformação logarítmica comprime diferenças em genes de baixa expressão. Genes com contagem normalizada < 1 têm valores negativos após transformação.

6. **Gene-level apenas:** Os dados são de expressão gênica total (*gene-level*), não de isoformas. Inferências sobre isoformas específicas (ex.: PKM1 vs. PKM2) não são suportadas.

### Limitações das comparações entre coortes

7. **Confundimento condição-coorte (GTEx):** Todos os tumores KIRP provêm do TCGA; todos os controles GTEx são de autópsias. O efeito tumoral e o efeito de coorte (diferenças de coleta, processamento, idade, causa de morte, composição celular, isquemia, profundidade de sequenciamento) são perfeitamente confundidos e não separáveis. A análise com GTEx é **estritamente exploratória** e seus resultados devem ser interpretados com extrema cautela.

8. **Controles não independentes:** As duas análises com tecido adjacente (pareada e não-pareada KIRP) compartilham os mesmos 32 controles. A concordância entre elas (CCC = 0,974) **não constitui validação independente**. O tecido adjacente pode apresentar alterações de campo, inflamação ou contaminação tumoral e não é necessariamente tecido renal completamente normal.

### Limitações da generalização

9. **Heterogeneidade molecular do KIRP:** O KIRP é molecularmente heterogêneo. A classificação OMS 2022 não utiliza mais a dicotomia tipo 1/tipo 2; entidades como carcinoma renal deficiente em fumarato hidratase (FH) são categorias distintas. A análise agrupada (*pooled*) não captura diferenças entre subtipos moleculares. Os resultados referem-se à média dos 288 tumores KIRP disponíveis e não devem ser generalizados para subtipos específicos.

10. **Escopo restrito:** A concordância entre comparadores foi avaliada apenas para os 106 genes do metabolismo central do carbono. A generalização para outros conjuntos gênicos ou para o transcriptoma completo requer verificação independente para cada contexto.

11. **Generalização para outros tipos de câncer renal:** Este estudo analisou exclusivamente KIRP. Os achados não se aplicam a carcinomas de células claras (KIRC), cromófobos (KICH) ou outros subtipos.

### Limitações computacionais

12. **Sensibilidade a parâmetros:** O resultado do teste camera para hsa00020 (FDR = 0,0012 com correlação predefinida de 0,01) não se manteve na análise de sensibilidade com correlação estimada por conjunto (FDR = 0,437 com `inter.gene.cor = NA`). A escolha do parâmetro de correlação intergênica afeta substancialmente as conclusões.

13. **Dependência de anotações externas:** As anotações KEGG (release 119.0, julho/2026) foram congeladas em arquivo versionado. Atualizações futuras do KEGG podem alterar a composição dos conjuntos gênicos.

14. **Reprodutibilidade dos HTMLs 3D:** As visualizações 3D dependem da biblioteca plotly.js (versão 2.25.2). Versões futuras podem alterar o comportamento de renderização. Os arquivos `*_files/` devem ser mantidos junto aos HTMLs.

---

## Declaração de Uso de Inteligência Artificial

Ferramentas de inteligência artificial generativa foram utilizadas como recursos de apoio ao longo do desenvolvimento deste projeto. Em nenhum momento as ferramentas substituíram o julgamento científico dos autores, a interpretação dos resultados ou a responsabilidade pela integridade acadêmica do trabalho. Todas as decisões metodológicas, análises estatísticas, interpretações biológicas e conclusões são de responsabilidade exclusiva dos autores.

### Sumário de uso de IA por etapa

| Ferramenta de IA | Fabricante | Etapas de uso | Natureza da contribuição |
|:---|:---|:---|:---|
| **ChatGPT-5.5** | OpenAI | Redação e revisão do manuscrito; estruturação argumentativa; formatação ABNT/Vancouver; tradução | Assistência na redação científica, organização de seções, verificação de consistência textual e referencial. Nenhum dado, análise ou resultado foi gerado por esta ferramenta. |
| **Codex (OpenAI)** | OpenAI | Desenvolvimento e depuração de scripts R e Python; implementação de funções do limma, ggplot2, plotly, igraph; criação de visualizações 3D; auditoria de código | Geração de código sob supervisão humana. Todo código foi revisado, testado e validado pelos autores antes da execução. O pipeline foi executado integralmente em ambiente controlado com sessão registrada. |
| **DeepSeek-v4-pro** | Hangzhou DeepSeek Artificial Intelligence | Verificação de consistência referencial; revisão cruzada de valores numéricos entre scripts, resultados e manuscrito; sugestões de melhorias metodológicas | Auditoria de consistência entre manuscrito e outputs. Identificação de incongruências numéricas e sugestões de correção. Todas as alterações foram validadas pelos autores. |
| **GitHub Copilot** | GitHub/Microsoft | Sugestões pontuais de código durante edição no IDE | Complementação de código boilerplate. Nenhum bloco completo de análise foi gerado exclusivamente por esta ferramenta. |
| **Claude (Anthropic)** | Anthropic | Revisão de estilo e clareza do manuscrito; verificação de conformidade com diretrizes de periódicos | Revisão de linguagem e aderência a normas de publicação. |

### Transparência adicional

- **Nenhum dado foi gerado sinteticamente:** Todos os dados são públicos do TCGA/GTEx via UCSC Xena.
- **Nenhuma análise estatística foi delegada à IA:** Modelos lineares (limma), testes de enriquecimento e métricas de concordância foram implementados em scripts R versionados, revisados e executados pelos autores.
- **Sessão registrada:** O ambiente computacional completo está documentado em `results/v3/sessionInfo.txt` e `environment/packages.csv`.
- **Código aberto e auditável:** Todos os scripts estão disponíveis neste repositório sob licença MIT.

---

## Como Citar

### Formato ABNT (NBR 6023:2018)

LOUREIRO, Kamila da Conceição; SANTOS, Ryan de Paulo; FREITAS, Letícia Maria Dias; SILVA, Ivine Souza; PECLY, Maria Eduarda Peixoto Soares. **KIRP-Glycolysis-Transcriptomics**: perfil transcriptômico do metabolismo central do carbono no carcinoma renal papilar. Versão 3.1.0. [S. l.], 2026. Código-fonte. Disponível em: https://github.com/santosry/KIRP-Glycolysis-Transcriptomics. Acesso em: [data de acesso].

### Formato sugerido (software)

Santos, R. P., Loureiro, K. C., Freitas, L. M. D., Silva, I. S., & Pecly, M. E. P. S. (2026). *KIRP-Glycolysis-Transcriptomics* (Version 3.1.0) [Computer software]. https://github.com/santosry/KIRP-Glycolysis-Transcriptomics

### Metadados

Consulte o arquivo [CITATION.cff](CITATION.cff) para metadados de citação no formato Citation File Format (CFF 1.2.0), compatível com GitHub, Zenodo e Zotero.

---

## Referências

As referências completas no formato ABNT (NBR 6023:2018) com DOI, link de acesso e data de acesso:

1. LINEHAN, W. M. et al. Comprehensive molecular characterization of papillary renal-cell carcinoma. **New England Journal of Medicine**, v. 374, n. 2, p. 135-145, 2016. DOI: [10.1056/NEJMoa1505917](https://doi.org/10.1056/NEJMoa1505917). Disponível em: https://www.nejm.org/doi/full/10.1056/NEJMoa1505917. Acesso em: 11 jul. 2026.

2. MOCH, H. et al. The 2022 World Health Organization classification of tumours of the urinary system and male genital organs — part A: renal, penile, and testicular tumours. **European Urology**, v. 82, n. 5, p. 469-482, 2022. DOI: [10.1016/j.eururo.2022.07.011](https://doi.org/10.1016/j.eururo.2022.07.011). Disponível em: https://www.sciencedirect.com/science/article/pii/S0302283822025074. Acesso em: 11 jul. 2026.

3. VANDER HEIDEN, M. G.; CANTLEY, L. C.; THOMPSON, C. B. Understanding the Warburg effect: the metabolic requirements of cell proliferation. **Science**, v. 324, n. 5930, p. 1029-1033, 2009. DOI: [10.1126/science.1160809](https://doi.org/10.1126/science.1160809). Disponível em: https://www.science.org/doi/10.1126/science.1160809. Acesso em: 11 jul. 2026.

4. PATRA, K. C.; HAY, N. The pentose phosphate pathway and cancer. **Trends in Biochemical Sciences**, v. 39, n. 8, p. 347-354, 2014. DOI: [10.1016/j.tibs.2014.06.005](https://doi.org/10.1016/j.tibs.2014.06.005). Disponível em: https://www.cell.com/trends/biochemical-sciences/fulltext/S0968-0004(14)00097-2. Acesso em: 11 jul. 2026.

5. TESLAA, T. et al. The pentose phosphate pathway in health and disease. **Nature Metabolism**, v. 5, n. 8, p. 1275-1289, 2023. DOI: [10.1038/s42255-023-00863-2](https://doi.org/10.1038/s42255-023-00863-2). Disponível em: https://www.nature.com/articles/s42255-023-00863-2. Acesso em: 11 jul. 2026.

6. MARTINEZ-REYES, I.; CHANDEL, N. S. Coupling Krebs cycle metabolites to signalling in immunity and cancer. **Nature Metabolism**, v. 1, n. 1, p. 16-33, 2019. DOI: [10.1038/s42255-018-0014-7](https://doi.org/10.1038/s42255-018-0014-7). Disponível em: https://www.nature.com/articles/s42255-018-0014-7. Acesso em: 11 jul. 2026.

7. SANCHEZ, D. J.; SIMON, M. C. Metabolic alterations in hereditary and sporadic renal cell carcinoma. **Nature Reviews Nephrology**, v. 20, n. 8, p. 521-534, 2024. DOI: [10.1038/s41581-024-00821-3](https://doi.org/10.1038/s41581-024-00821-3). Disponível em: https://www.nature.com/articles/s41581-024-00821-3. Acesso em: 11 jul. 2026.

8. GE, X. et al. Fumarate inhibits PTEN to promote tumorigenesis and therapeutic resistance of type 2 papillary renal cell carcinoma. **Molecular Cell**, v. 82, n. 10, p. 1929-1944.e8, 2022. DOI: [10.1016/j.molcel.2022.03.011](https://doi.org/10.1016/j.molcel.2022.03.011). Disponível em: https://www.cell.com/molecular-cell/fulltext/S1097-2765(22)00235-0. Acesso em: 11 jul. 2026.

9. KOPPULA, P. et al. Cystine transporter regulation of pentose phosphate pathway dependency and disulfide stress exposes a targetable metabolic vulnerability in cancer. **Nature Cell Biology**, v. 22, n. 4, p. 476-486, 2020. DOI: [10.1038/s41556-020-0496-x](https://doi.org/10.1038/s41556-020-0496-x). Disponível em: https://www.nature.com/articles/s41556-020-0496-x. Acesso em: 11 jul. 2026.

10. GOLDMAN, M. J. et al. Visualizing and interpreting cancer genomics data via the Xena platform. **Nature Biotechnology**, v. 38, n. 6, p. 675-678, 2020. DOI: [10.1038/s41587-020-0546-8](https://doi.org/10.1038/s41587-020-0546-8). Disponível em: https://www.nature.com/articles/s41587-020-0546-8. Acesso em: 11 jul. 2026.

11. GTEX CONSORTIUM. The Genotype-Tissue Expression (GTEx) project. **Nature Genetics**, v. 45, n. 6, p. 580-585, 2013. DOI: [10.1038/ng.2653](https://doi.org/10.1038/ng.2653). Disponível em: https://www.nature.com/articles/ng.2653. Acesso em: 11 jul. 2026.

12. RITCHIE, M. E. et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. **Nucleic Acids Research**, v. 43, n. 7, p. e47, 2015. DOI: [10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007). Disponível em: https://academic.oup.com/nar/article/43/7/e47/2414268. Acesso em: 11 jul. 2026.

13. VIVIAN, J. et al. Toil enables reproducible, open source, big biomedical data analyses. **Nature Biotechnology**, v. 35, n. 4, p. 314-316, 2017. DOI: [10.1038/nbt.3772](https://doi.org/10.1038/nbt.3772). Disponível em: https://www.nature.com/articles/nbt.3772. Acesso em: 11 jul. 2026.

14. WU, D.; SMYTH, G. K. Camera: a competitive gene set test accounting for inter-gene correlation. **Nucleic Acids Research**, v. 40, n. 17, p. e133, 2012. DOI: [10.1093/nar/gks461](https://doi.org/10.1093/nar/gks461). Disponível em: https://academic.oup.com/nar/article/40/17/e133/2411347. Acesso em: 11 jul. 2026.

15. UHLÉN, M. et al. Tissue-based map of the human proteome. **Science**, v. 347, n. 6220, p. 1260419, 2015. DOI: [10.1126/science.1260419](https://doi.org/10.1126/science.1260419). Disponível em: https://www.science.org/doi/10.1126/science.1260419. Acesso em: 11 jul. 2026.

16. ARAN, D.; SIROTA, M.; BUTTE, A. J. Systematic pan-cancer analysis of tumour purity. **Nature Communications**, v. 6, p. 8971, 2015. DOI: [10.1038/ncomms9971](https://doi.org/10.1038/ncomms9971). Disponível em: https://www.nature.com/articles/ncomms9971. Acesso em: 11 jul. 2026.

---

## Licença

Este projeto está licenciado sob a Licença MIT — veja o arquivo [LICENSE](LICENSE) para detalhes.

Copyright © 2026 Kamila da Conceição Loureiro, Ryan de Paulo Santos, Letícia Maria Dias Freitas, Ivine Souza Silva, Maria Eduarda Peixoto Soares Pecly.

---

*Última atualização: 12 de julho de 2026. Versão 3.1.0.*
