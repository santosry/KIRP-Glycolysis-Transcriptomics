# KIRP Glycolysis · PPP · TCA — Transcriptomics

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R 4.6.0](https://img.shields.io/badge/R-4.6.0-276DC3)](https://www.r-project.org/)
[![v3.1.0](https://img.shields.io/badge/version-3.1.0-green)]()
[![Python 3](https://img.shields.io/badge/Python-3.x-3776AB)](https://www.python.org/)

**Transcriptoma das vias glicolítica, das pentoses fosfato e do ciclo do ácido cítrico no carcinoma renal papilar: análise pareada e sensibilidade ao tecido de referência.**

> ⚠️ **Nota sobre a abrangência das vias KEGG:** A via hsa00010 (Glicólise/Gliconeogênese) inclui, na anotação KEGG *release* 119.0, famílias gênicas cuja função primária é o metabolismo de etanol e a detoxificação de aldeídos (ADH: *ADH1A/B/C, ADH4/5/6/7*; ALDH: *ALDH1B1, ALDH2, ALDH3A1/A2/B1/B2, ALDH7A1, ALDH9A1*). Estes ~25 genes são co-anotados na via porque o etanol pode ser convertido a acetil-CoA, mas **não pertencem ao metabolismo central do carbono em sentido estrito**. Sua inclusão decorre da estratégia de usar a anotação KEGG completa sem curadoria manual — decisão que privilegia reprodutibilidade sobre pureza funcional. A discussão científica concentra-se nos genes de função estabelecida nas três vias (HK2, ALDOB, FBP1, PCK1, G6PD, TKT, ENO2 e demais enzimas *core*). Os resultados de ADH/ALDH constam integralmente nos dados suplementares (Supplementary Table S1).

---

## Índice

1. [Questão Científica](#questão-científica)
2. [Vias Analisadas](#vias-analisadas)
3. [Por que estas três vias?](#por-que-estas-três-vias)
4. [Principais Achados](#principais-achados)
5. [Desenho dos Comparadores](#desenho-dos-comparadores)
6. [Dados](#dados)
7. [Estrutura do Repositório](#estrutura-do-repositório)
8. [Reprodução](#reprodução)
9. [Instruções para Agentes de IA](#instruções-para-agentes-de-ia)
10. [Visualizações 3D Interativas](#visualizações-3d-interativas)
11. [Limitações](#limitações)
12. [Declaração de Uso de Inteligência Artificial](#declaração-de-uso-de-inteligência-artificial)
13. [Como Citar](#como-citar)
14. [Referências](#referências)
15. [Licença](#licença)

---

## Questão Científica

Qual é o perfil transcriptômico dos genes das vias **glicolítica (hsa00010), das pentoses fosfato (hsa00030) e do ciclo do ácido cítrico (hsa00020)** no carcinoma renal papilar (KIRP) quando analisado com tecido normal adjacente pareado como comparador primário, e qual a sensibilidade dos resultados à escolha do tecido de referência?

---

## Vias Analisadas

O estudo investiga **três vias do metabolismo central do carbono**, tratadas com igual importância. Nenhuma via é hierarquicamente "principal": as três são analisadas em conjunto a partir de um modelo transcriptômico único sobre 31.633 genes, com extração *a posteriori* dos 106 genes de interesse — estratégia que evita viés de seleção e circularidade metodológica.

| KEGG ID | Via / Pathway | Genes no KEGG | Genes na Matriz |
|---------|---------------|:-------------:|:---------------:|
| hsa00010 | Glicólise / Gliconeogênese | 67 | 64 |
| hsa00030 | Via das Pentoses Fosfato | 31 | 30 |
| hsa00020 | Ciclo do Ácido Cítrico (TCA) | 30 | 29 |
| **União (genes únicos)** | **Metabolismo Central do Carbono** | **110** | **106** |

**Genes ausentes da matriz (n = 3):** *PDHA2*, *PGK2*, *PRPS1L1* — removidos durante a filtragem de baixa expressão (> 90% das amostras com expressão ≤ 1). Análise de sensibilidade confirmou que nenhum gene das três vias foi afetado exclusivamente por este critério. Matriz final: **106 genes**.

> ⚠️ **Alcance da anotação KEGG hsa00010:** Aproximadamente 25 dos 64 genes desta via são álcool desidrogenases (ADH) e aldeído desidrogenases (ALDH), enzimas do metabolismo de etanol e detoxificação de aldeídos, co-anotadas porque o etanol pode entrar como acetil-CoA. Estes genes **não pertencem ao metabolismo central do carbono em sentido estrito**. São mantidos na análise por fidelidade à anotação KEGG *release* 119.0 (sem curadoria manual — decisão que privilegia reprodutibilidade). A discussão concentra-se nas enzimas *core* (HK2, ALDOB, FBP1, PCK1, G6PD, TKT, ENO2 e demais). Resultados completos de ADH/ALDH em Supplementary Table S1.

---

## Por que estas três vias?

### Glicólise / Gliconeogênese (hsa00010)

A via glicolítica é o ponto de partida do metabolismo central do carbono: converte glicose em piruvato, gerando ATP e intermediários biossintéticos que alimentam todas as demais vias. A gliconeogênese, via reversa, é crítica no rim: o túbulo proximal renal é um dos principais sítios gliconeogênicos do organismo humano, expressando *FBP1*, *PCK1* e *ALDOB* em níveis basais elevados. No câncer, o efeito Warburg — aumento da captação de glicose e fermentação a lactato mesmo na presença de oxigênio — é uma das alterações metabólicas mais documentadas (VANDER HEIDEN; CANTLEY; THOMPSON, 2009). A inclusão desta via no KIRP justifica-se por duas razões: (i) a perda de enzimas gliconeogênicas no tumor serve como marcador de desdiferenciação epitelial tubular; (ii) isoenzimas como *HK2* e *PKM2* são alvos terapêuticos em múltiplos cânceres e sua expressão no KIRP não havia sido caracterizada em análise pareada.

### Via das Pentoses Fosfato (hsa00030)

A PPP ramifica-se da glicólise na glicose-6-fosfato e desempenha duas funções essenciais para a proliferação celular: (i) produção de **NADPH**, o principal redutor celular que sustenta a defesa antioxidante (glutationa, tiorredoxina) e a biossíntese de lipídeos; (ii) produção de **ribose-5-fosfato**, precursor da síntese de nucleotídeos (DNA, RNA). Em carcinomas renais, a PPP é particularmente relevante porque a enzima limitante **G6PD** é regulada pela via NRF2 — fator de transcrição frequentemente ativado em KIRP com deficiência de FH (KOPPULA *et al.*, 2020; PATRA; HAY, 2014). A inclusão desta via permite testar se há evidência transcriptômica de ativação coordenada do braço oxidativo da PPP.

### Ciclo do Ácido Cítrico — TCA (hsa00020)

O ciclo de Krebs é a encruzilhada central do metabolismo oxidativo: oxida acetil-CoA a CO₂, gera NADH e FADH₂ para a cadeia respiratória e fornece intermediários para biossíntese de aminoácidos, lipídeos e heme. No KIRP, o TCA tem relevância única: mutações no gene ***FH*** (fumarato hidratase) são a base molecular do subtipo hereditário com deficiência de FH (HLRCC), e o acúmulo de fumarato promove succinação de PTEN conectando diretamente o ciclo à tumorigênese (GE *et al.*, 2022). Adicionalmente, genes como *SDHA-D*, *IDH1/2* e *MDH2* codificam enzimas mitocondriais cujas mutações são recorrentes em neoplasias renais (SANCHEZ; SIMON, 2024). A inclusão sistemática dos 29 genes do TCA permite testar formalmente se o ciclo como conjunto gênico apresenta alteração coordenada — hipótese avaliada via *camera*.

### Integração das três vias

As três vias não operam isoladamente: a glicólise alimenta a PPP (via glicose-6-fosfato) e o TCA (via piruvato → acetil-CoA); a PPP fornece NADPH para biossíntese redutora; o TCA fornece citrato para lipogênese e intermediários para aminoácidos. Analisar as três vias em conjunto — 106 genes extraídos *a posteriori* de um modelo transcriptômico único — evita o viés de seleção que ocorreria ao restringir previamente o universo gênico apenas a essas vias. A análise conjunta também revela **genes compartilhados** entre vias (7 dos 35 DEGs pertencem a duas ou três vias), evidenciando pontos de conexão metabólica potencialmente relevantes.

---

## Principais Achados (v3.1.0)

### Análise pareada primária (32 pares tumor-adjacente KIRP)

- **35/106 genes (33,0%)** diferencialmente expressos (|log₂FC| > 1, FDR < 0,05): 9 aumentados, 26 diminuídos
- 7 genes compartilhados entre duas ou três vias metabólicas (Tabela abaixo)

#### Genes compartilhados entre vias (DEGs de dupla anotação KEGG)

| Gene | Vias | log₂FC | FDR | Dir. | Função |
|------|------|:------:|:---:|:----:|--------|
| *ALDOA* | hsa00010 + hsa00030 | +1,12 | 1,0×10⁻⁶ | ▲ | Aldolase A (glicólise) |
| *ALDOB* | hsa00010 + hsa00030 | −8,66 | 4,3×10⁻¹³ | ▼ | Aldolase B (gliconeogênese renal) |
| *ALDOC* | hsa00010 + hsa00030 | −1,08 | 1,5×10⁻² | ▼ | Aldolase C (cerebral) |
| *FBP1* | hsa00010 + hsa00030 | −3,07 | 1,1×10⁻⁸ | ▼ | Frutose-1,6-bisfosfatase 1 |
| *FBP2* | hsa00010 + hsa00030 | −1,57 | 4,7×10⁻⁸ | ▼ | Frutose-1,6-bisfosfatase 2 |
| *PCK1* | hsa00010 + hsa00020 | −5,11 | 1,7×10⁻⁹ | ▼ | PEP carboxiquinase 1 |
| *PCK2* | hsa00010 + hsa00020 | −1,69 | 6,2×10⁻⁵ | ▼ | PEP carboxiquinase 2 |

Dos 7 genes: 5 pertencem simultaneamente à glicólise e à PPP (ALDOA/B/C, FBP1/2 — enzimas que atuam no ponto de ramificação frutose-1,6-bisfosfato); 2 pertencem à glicólise e ao TCA (PCK1/2 — enzimas que conectam oxaloacetato à gliconeogênese).

| Via | Genes | DEGs | ▲ Up | ▼ Down |
|-----|:-----:|:----:|:----:|:------:|
| Glicólise / Gliconeogênese (hsa00010) | 64 | 26 | 7 | 19 |
| Via das Pentoses Fosfato (hsa00030) | 30 | 12 | 3 | 9 |
| Ciclo do Ácido Cítrico (hsa00020) | 29 | 4 | 0 | 4 |
| **Total (genes únicos)** | **106** | **35** | **9** | **26** |

### Genes de maior magnitude

| Gene | Via(s) | log₂FC | IC 95% | FDR | Direção |
|------|--------|:------:|--------|:---:|:-------:|
| *ALDOB* | hsa00010, hsa00030 | −8,66 | [−9,96; −7,36] | 4,3×10⁻¹³ | ↓ |
| *ADH1C* | hsa00010 | −6,01 | [−7,06; −4,95] | 1,9×10⁻¹⁵ | ↓ |
| *ALDH3B2* | hsa00010 | −5,19 | [−5,93; −4,46] | 8,5×10⁻²¹ | ↓ |
| *PCK1* | hsa00010, hsa00020 | −5,11 | [−6,53; −3,69] | 1,7×10⁻⁹ | ↓ |
| *HK2* | hsa00010 | +3,34 | [2,51; 4,17] | 1,9×10⁻¹² | ↑ |
| *FBP1* | hsa00010, hsa00030 | −3,07 | [−3,90; −2,23] | 1,1×10⁻⁸ | ↓ |
| *TKT* | hsa00030 | +1,60 | [1,19; 2,00] | 8,1×10⁻¹² | ↑ |
| *G6PD* | hsa00030 | +1,49 | [0,89; 2,09] | 1,8×10⁻⁶ | ↑ |

### Teste de conjuntos gênicos (camera)

| Via | NGenes | Direção | P | FDR |
|-----|:------:|:-------:|:---:|:---:|
| hsa00020 (TCA) | 29 | ↓ Down | 0,00039 | **0,159 (não significativo)** |
| hsa00010 (Glicólise) | 64 | ↓ Down | 0,968 | 0,968 |
| hsa00030 (PPP) | 30 | ↑ Up | 0,248 | 0,352 |

Nenhuma das três vias apresentou alteração coordenada estatisticamente significativa após correção para testes múltiplos. O TCA apresentou o menor valor nominal (FDR = 0,159).

### Concordância entre comparadores

| Métrica | Pareado vs. TCGA-KIRP expandida | Pareado vs. GTEx |
|---------|:-------------------------------:|:----------------:|
| CCC de Lin | 0,990 | 0,795 |
| MAE (log₂FC) | 0,19 | 0,94 |
| Viés (pareado − comparador) | +0,03 | −0,79 |
| Limites de concordância 95% | [−0,34; 0,40] | [−2,24; 0,67] |
| Concordância direcional | 100% (75/75) | 60,2% (59/98) |
| Genes discordantes (todos os 106) | 4 | 45 |

> ⚠️ **GTEx como referência externa introduz viés sistemático** de −0,79 e reverte a direção em 45 genes. As duas análises com tecido adjacente KIRP compartilham os mesmos 32 controles; a alta concordância **não** constitui validação independente.

---

## Desenho dos Comparadores

| Comparador | Delineamento | N | Papel |
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
| **Amostras** | 445 (288 KIRP, 28 GTEx, 32 adjacentes KIRP, 97 outros TCGA) |
| **Genes (brutos)** | 58.581 |
| **Genes (após filtragem)** | 31.633 |
| **Escala** | log₂(norm_count + 1) — pré-transformada pelo UCSC Xena |
| **Faixa** | 0 a 20,1 (mediana = 11,15) |
| **Acesso** | 11 de julho de 2026 |

### Como obter os dados

#### Arquivo kidney.tsv (~220 KB, versionado)

Este arquivo está em `data/raw/kidney.tsv`. Para baixá-lo diretamente do UCSC Xena:

1. Acesse [https://xenabrowser.net/datapages/](https://xenabrowser.net/datapages/)
2. Selecione: **GDC TCGA Kidney Papillary Cell Carcinoma (KIRP)** + **GTEx Kidney**
3. Escolha: **gene expression RNAseq — RSEM norm_count (log₂)**
4. Salve como `data/raw/kidney.tsv`

#### Transcriptoma completo (~178 MB, NÃO versionado)

```bash
python3 scripts/download_full_transcriptome.py
```

Verificação:
```bash
sha256sum data/raw/kidney.tsv
# Esperado: 29154696504ce365e681d9c319fe352a6c84c5ae87798c8f1cce3c11159f7ea2
```

---

## Estrutura do Repositório

```
KIRP-Glycolysis-Transcriptomics/
│
├── data/
│   ├── raw/
│   │   ├── kidney.tsv                    # Painel amostral (~220 KB, versionado)
│   │   └── kidney_transcriptome.tsv      # Transcriptoma completo (~178 MB, NÃO versionado)
│   ├── processed/                        # Matrizes processadas (RDS)
│   │   ├── expression_hsa00010.rds       # Expressão — Glicólise
│   │   ├── expression_hsa00020.rds       # Expressão — TCA
│   │   ├── expression_hsa00030.rds       # Expressão — PPP
│   │   └── expression_matrix_full.rds    # Matriz completa
│   ├── metadata/                         # Metadados de genes
│   └── provenance/                       # Registro de proveniência
│
├── run_all.R                             # Executa pipeline completo (01→06)
├── scripts/
│   ├── 00_environment.R                  # Verificação de dependências
│   ├── 01_load_data.R                    # Carrega kidney.tsv, filtra, metadados
│   ├── 02_differential_expression.R      # limma: paired + TCGA + GTEx
│   ├── 03_volcano_plots.R                # Volcano PNG + 3D HTML por via
│   ├── 04_ppi_network.R                  # Rede coexpressão (PNG + 3D HTML)
│   ├── 05_functional_enrichment.R        # camera + ORA KEGG (Up/Down)
│   └── 06_concordance_and_tables.R       # CCC, Bland-Altman, tabelas S1/DEG
│
├── results/
│   └── v3/                               # Todos os outputs da versão 3
│       ├── tables/
│       │   ├── Supplementary_Table_S1.csv    # Tabela suplementar (106 genes)
│       │   ├── DEG_hsa00010.csv              # DEGs — Glicólise
│       │   ├── DEG_hsa00020.csv              # DEGs — TCA
│       │   ├── DEG_hsa00030.csv              # DEGs — PPP
│       │   ├── discordant_genes.csv          # Discordantes (Paired vs GTEx)
│       │   ├── camera_gene_sets.csv          # Resultados camera
│       │   ├── comparator_concordance.csv    # Métricas de concordância
│       │   ├── ppi_3d_centrality.csv         # Centralidade da rede 3D
│       │   └── ppi_3d_summary.csv            # Resumo da rede 3D
│       ├── figures/
│       │   ├── Volcano_hsa00010.png / _3D.html   # Volcano — Glicólise
│       │   ├── Volcano_hsa00020.png / _3D.html   # Volcano — TCA
│       │   ├── Volcano_hsa00030.png / _3D.html   # Volcano — PPP
│       │   ├── PPI_network_3D.html                # Rede 3D interativa
│       │   ├── BlandAltman_Paired_vs_GTEx.png     # Bland-Altman
│       │   ├── PCA_transcriptome.png              # PCA transcriptoma
│       │   ├── plotSA_paired.png                  # Diagnóstico limma
│       │   ├── QC_density.png                     # Densidade de expressão
│       │   └── Paired_*.png                       # Gráficos pareados por gene
│       ├── supplementary/                  # Material suplementar (S2-S8)
│       ├── audit/                          # Relatórios de auditoria
│       ├── sessionInfo.txt                 # Sessão R
│       ├── checksums_sha256.txt            # Checksums SHA256
│       └── results_manifest.csv            # Manifesto de outputs
│
├── tests/
│   └── testthat/test_pipeline.R            # Testes automatizados
│
├── environment/
│   ├── packages.csv                        # Versões de pacotes R
│   └── sessionInfo.txt                     # Sessão R completa
│
├── renv.lock                               # Lockfile do renv
├── .github/workflows/ci.yml                # CI (GitHub Actions)
├── .gitignore
├── .Rprofile
├── LICENSE                                 # MIT
├── CITATION.cff                            # Metadados de citação
├── VERSION                                 # v3.1.0
├── README.md                               # Este arquivo
├── REPRODUCIBILITY.md                      # Guia de reprodutibilidade
└── DATA_PROVENANCE.md                      # Proveniência dos dados
```

---

## Reprodução

### Pré-requisitos

- **R** ≥ 4.6.0 (Bioconductor 3.23)
- **Python** ≥ 3.8 (apenas para download dos dados)
- **Pacotes R:** `limma`, `dplyr`, `rio`, `ggplot2`, `ggrepel`, `igraph`, `plotly`, `htmlwidgets`, `clusterProfiler`, `org.Hs.eg.db` e dependências (`environment/packages.csv`)
- **Conexão com internet:** necessária para download dos dados (~2 GB durante o processo)
- **Espaço em disco:** ~2 GB (download + processamento); ~250 MB (resultados finais)
- **Memória RAM:** ≥ 8 GB

### Pipeline completo

```bash
# 1. Clonar
git clone https://github.com/santosry/KIRP-Glycolysis-Transcriptomics.git
cd KIRP-Glycolysis-Transcriptomics

# 2. O arquivo data/raw/kidney.tsv já está versionado (110 genes, 445 amostras)
#    Para baixar do UCSC Xena, veja instruções na seção 'Dados'

# 3. Executar pipeline completo
Rscript run_all.R

# 4. Testes
R -e "testthat::test_dir('tests/testthat')"
```

### Verificação de integridade

```bash
cd results/v3
sha256sum -c checksums_sha256.txt
```

---

## Instruções para Agentes de IA

> **Esta seção é destinada a agentes de IA (Codex, Claude, ChatGPT, DeepSeek, Copilot, etc.) que executem este repositório. Leia antes de qualquer comando.**

### Regras gerais

1. **Diretório de trabalho:** Executar scripts a partir da raiz do repositório (`KIRP-Glycolysis-Transcriptomics/`).

2. **RENV:** Verificar sincronia antes de scripts R:
   ```r
   renv::status()
   renv::restore()  # se necessário
   ```
   Se renv indisponível, instalar pacotes listados em `environment/packages.csv`. **Não modificar `renv.lock`** sem solicitação explícita.

3. **Ordem de execução obrigatória:**
   ```
   1. scripts/00_environment.R                 (verifica dependências)
   2. run_all.R                                (pipeline completo 01→06)
   ```
   Ou passo a passo:
   ```
   1. scripts/01_load_data.R
   2. scripts/02_differential_expression.R
   3. scripts/03_volcano_plots.R
   4. scripts/04_ppi_network.R
   5. scripts/05_functional_enrichment.R
   6. scripts/06_concordance_and_tables.R
   ```

4. **Dados obrigatórios:**
   - `data/raw/kidney.tsv` — **versionado** (~220 KB). Contém IDs das 445 amostras renais. **Não sobrescrever**.
   - `data/raw/kidney_transcriptome.tsv` — **NÃO versionado** (~178 MB). Baixar via `download_full_transcriptome.py`.

5. **Timeout:** Download do transcriptoma completo leva 10-20 min. Configurar timeout ≥ 1200 s.

6. **Não executar** scripts com prefixo `10_` (STRING) sem solicitação explícita — dependem de API externa e ~500 MB de download. A rede PPI 3D atual (`16b_ppi_3d_correlation.R`) usa correlação de Pearson (offline, 100% reprodutível).

7. **Scripts v2 (prefixos 01-15):** legado da versão 2 (apenas glicólise, 66 genes). Usar `pipeline_v3.R` e `16_*.R` para a versão atual (3 vias, 106 genes).

8. **Não modificar** arquivos em `results/v3/` manualmente.

9. **Manuscrito** (`.Rmd`, `.pdf`): **não versionar**, não publicar.

10. **Ao reportar resultados:** sempre citar números exatos dos arquivos de saída (`results/v3/tables/`), não do README.

11. **Após gerar resultados:** verificar checksums (`sha256sum -c results/v3/checksums_sha256.txt`) e executar testes (`testthat::test_dir("tests/testthat")`).

### Checklist de verificação

- [ ] `data/raw/kidney.tsv` existe (~220 KB)
- [ ] `data/raw/kidney_transcriptome.tsv` existe (~178 MB) ou foi baixado com sucesso
- [ ] `renv::status()` retorna sincronizado
- [ ] `source("scripts/00_environment.R")` sem erros
- [ ] `pipeline_v3.R` executa até o fim
- [ ] `pipeline_v3_addendum.R` executa até o fim
- [ ] `results/v3/tables/Supplementary_Table_S1.csv` existe
- [ ] `results/v3/figures/Volcano_*_3D.html` e `PPI_network_3D.html` existem
- [ ] `testthat::test_dir("tests/testthat")` — todos passam
- [ ] Nenhum arquivo de manuscrito rastreado pelo Git

---

## Visualizações 3D Interativas

Abra os arquivos HTML em `results/v3/figures/` em qualquer navegador moderno:

### Volcano Plots 3D

| Arquivo | Via | Genes |
|---------|-----|:-----:|
| `Volcano_hsa00010_3D.html` | Glicólise / Gliconeogênese | 64 |
| `Volcano_hsa00020_3D.html` | Ciclo do Ácido Cítrico (TCA) | 29 |
| `Volcano_hsa00030_3D.html` | Via das Pentoses Fosfato | 30 |

**Eixos:** X = log₂FC | Y = −log₁₀(FDR) | Z = AveExpr

**Cores:** 🔵 Azul = aumentado | 🟣 Roxo = diminuído | ⚪ Cinza = NS

**Interações:** arrastar (rotacionar) | scroll (zoom) | duplo-clique (reset) | hover (detalhes)

### Rede de Coexpressão 3D

| Arquivo | Método | Threshold |
|---------|--------|:---------:|
| `PPI_network_3D.html` | Pearson | \|r\| > 0,6 |

- **85 nós** no componente gigante, **935 arestas**
- Nós proporcionais ao degree | Cores por regulação
- Layout Fruchterman-Reingold 3D (seed = 42, determinístico)

---

## Limitações

### Desenho experimental

1. **32 pares apenas:** pequenos efeitos (|log₂FC| < 1) podem não ser detectados.
2. **Bulk RNA-seq:** não distingue tipos celulares. Redução de *ALDOB*, *FBP1*, *PCK1* pode refletir perda de células tubulares, não reprogramação metabólica. Sem deconvolução celular, esta hipótese concorrente não pode ser excluída.
3. **Sem validação externa:** resultados específicos da coorte TCGA-KIRP.

### Medição

4. **mRNA ≠ proteína ≠ atividade enzimática ≠ fluxo:** inferências sobre atividade metabólica são indiretas.
5. **Dados pré-transformados:** escala log₂(norm_count + 1) fornecida pelo Xena.
6. **Gene-level apenas:** inferências sobre isoformas (PKM1 vs PKM2) não suportadas.

### Comparações entre coortes

7. **Confundimento GTEx:** todos os tumores TCGA, todos os controles GTEx de autópsias. Efeito tumoral e de coorte são indissociáveis.
8. **Controles não independentes:** análises pareada e TCGA-KIRP expandida compartilham os mesmos 32 controles. Concordância não é validação independente.

### Generalização

9. **Heterogeneidade do KIRP:** OMS 2022 não usa mais dicotomia tipo 1/2. Análise agrupada não captura diferenças entre subtipos.
10. **Escopo restrito:** concordância avaliada apenas nos 106 genes destas três vias.
11. **Não se aplica** a KIRC, KICH ou outros subtipos renais.

### Computacionais

12. **Sensibilidade a parâmetros:** resultado do camera para TCA (FDR = 0,159 (não significativo)) não se manteve com *inter.gene.cor = NA* (FDR = 0,437).
13. **Anotações KEGG:** congeladas na release 119.0 (julho/2026). Atualizações futuras podem alterar a composição dos conjuntos gênicos.
14. **Genes periféricos em hsa00010:** ~25 dos 64 genes da via glicolítica são ADH/ALDH (metabolismo de etanol/detoxificação), co-anotados na via KEGG mas funcionalmente não pertencentes ao metabolismo central do carbono. Sua inclusão decorre da estratégia de usar a anotação KEGG completa sem curadoria manual. Resultados relativos a estes genes devem ser interpretados com cautela quanto à relevância para o metabolismo energético tumoral.

---

## Declaração de Uso de Inteligência Artificial

Ferramentas de IA generativa foram utilizadas como recursos de apoio, sem substituir o julgamento científico dos autores. Todas as decisões metodológicas, análises estatísticas, interpretações biológicas e conclusões são de responsabilidade exclusiva dos autores.

| Ferramenta | Fabricante | Etapas | Contribuição |
|:-----------|:-----------|:-------|:-------------|
| **ChatGPT-5.5** | OpenAI | Redação e revisão do manuscrito; estrutura argumentativa; ABNT/Vancouver | Assistência na redação científica. Nenhum dado ou análise foi gerado por esta ferramenta. |
| **Codex** | OpenAI | Desenvolvimento e depuração de scripts R/Python; funções limma, ggplot2, plotly, igraph; visualizações 3D; auditoria | Geração de código sob supervisão. Todo código revisado, testado e validado pelos autores. |
| **DeepSeek-v4-pro** | Hangzhou DeepSeek AI | Verificação cruzada de valores numéricos; auditoria de consistência manuscrito vs. outputs | Identificação de incongruências. Todas as correções validadas manualmente. |
- ✅ Nenhum dado gerado sinteticamente
- ✅ Nenhuma análise estatística delegada integralmente à IA
- ✅ Sessão registrada: `results/v3/sessionInfo.txt`, `environment/packages.csv`, `renv.lock`

---

## Como Citar

### ABNT (NBR 6023:2018)

LOUREIRO, Kamila da Conceição; SANTOS, Ryan de Paulo; FREITAS, Letícia Maria Dias; SILVA, Ivine Souza; PECLY, Maria Eduarda Peixoto Soares. **KIRP-Glycolysis-Transcriptomics**: perfil transcriptômico dos genes das vias glicolítica, das pentoses fosfato e do ciclo do ácido cítrico no carcinoma renal papilar (KIRP). Versão 3.1.0. [S. l.], 2026. Código-fonte. Disponível em: https://github.com/santosry/KIRP-Glycolysis-Transcriptomics. Acesso em: [data].

### Software

Santos, R. P., Loureiro, K. C., Freitas, L. M. D., Silva, I. S., & Pecly, M. E. P. S. (2026). *KIRP-Glycolysis-Transcriptomics: transcriptomic profiling of glycolysis, pentose phosphate pathway, and citrate cycle genes in papillary renal cell carcinoma* (Version 3.1.0) [Computer software]. https://github.com/santosry/KIRP-Glycolysis-Transcriptomics

Metadados completos em [CITATION.cff](CITATION.cff) (CFF 1.2.0, compatível com GitHub, Zenodo, Zotero).

---

## Referências

Referências em formato **ABNT (NBR 6023:2018)** com DOI, link de acesso e data de acesso:

1. CANCER GENOME ATLAS RESEARCH NETWORK. Comprehensive molecular characterization of papillary renal-cell carcinoma. **New England Journal of Medicine**, v. 374, n. 2, p. 135-145, 2016. DOI: [10.1056/NEJMoa1505917](https://doi.org/10.1056/NEJMoa1505917). Disponível em: https://www.nejm.org/doi/full/10.1056/NEJMoa1505917. Acesso em: 15 mar. 2026.

2. MOCH, H.; AMIN, M. B.; BERNEY, D. M.; COMPÉRAT, E. M.; GILL, A. J.; HARTMANN, A.; MENON, S.; RASPOLLINI, M. R.; RUBIN, M. A.; SRIGLEY, J. R.; TAN, P. H.; TICKOO, S. K.; TSUZUKI, T.; TURAJLIC, S.; CREE, I. A.; NETTO, G. J. The 2022 World Health Organization classification of tumours of the urinary system and male genital organs, part A: renal, penile, and testicular tumours. **European Urology**, v. 82, n. 5, p. 469-482, 2022. DOI: [10.1016/j.eururo.2022.07.011](https://doi.org/10.1016/j.eururo.2022.07.011). Disponível em: https://www.sciencedirect.com/science/article/pii/S0302283822025074. Acesso em: 28 mar. 2026.

3. VANDER HEIDEN, M. G.; CANTLEY, L. C.; THOMPSON, C. B. Understanding the Warburg effect: the metabolic requirements of cell proliferation. **Science**, v. 324, n. 5930, p. 1029-1033, 2009. DOI: [10.1126/science.1160809](https://doi.org/10.1126/science.1160809). Disponível em: https://www.science.org/doi/10.1126/science.1160809. Acesso em: 3 abr. 2026.

4. PATRA, K. C.; HAY, N. The pentose phosphate pathway and cancer. **Trends in Biochemical Sciences**, v. 39, n. 8, p. 347-354, 2014. DOI: [10.1016/j.tibs.2014.06.005](https://doi.org/10.1016/j.tibs.2014.06.005). Disponível em: https://www.cell.com/trends/biochemical-sciences/fulltext/S0968-0004(14)00097-2. Acesso em: 10 abr. 2026.

5. TESLAA, T.; RALSER, M.; FAN, J.; RABINOWITZ, J. D. The pentose phosphate pathway in health and disease. **Nature Metabolism**, v. 5, n. 8, p. 1275-1289, 2023. DOI: [10.1038/s42255-023-00863-2](https://doi.org/10.1038/s42255-023-00863-2). Disponível em: https://www.nature.com/articles/s42255-023-00863-2. Acesso em: 22 abr. 2026.

6. MARTINEZ-REYES, I.; CHANDEL, N. S. Coupling Krebs cycle metabolites to signalling in immunity and cancer. **Nature Metabolism**, v. 1, n. 1, p. 16-33, 2019. DOI: [10.1038/s42255-018-0014-7](https://doi.org/10.1038/s42255-018-0014-7). Disponível em: https://www.nature.com/articles/s42255-018-0014-7. Acesso em: 5 mai. 2026.

7. SANCHEZ, D. J.; SIMON, M. C. Metabolic alterations in hereditary and sporadic renal cell carcinoma. **Nature Reviews Nephrology**, v. 20, n. 8, p. 521-534, 2024. DOI: [10.1038/s41581-024-00821-3](https://doi.org/10.1038/s41581-024-00821-3). Disponível em: https://www.nature.com/articles/s41581-024-00821-3. Acesso em: 18 mai. 2026.

8. GE, X.; LYU, P.; GU, Y.; ZHOU, L.; LI, J.; XU, P.; LIU, J.; XU, Y.; LIU, Y.; WANG, Y.; GAN, B.; CHEN, J. Fumarate inhibits PTEN to promote tumorigenesis and therapeutic resistance of type 2 papillary renal cell carcinoma. **Molecular Cell**, v. 82, n. 10, p. 1929-1944.e8, 2022. DOI: [10.1016/j.molcel.2022.03.011](https://doi.org/10.1016/j.molcel.2022.03.011). Disponível em: https://www.cell.com/molecular-cell/fulltext/S1097-2765(22)00235-0. Acesso em: 2 jun. 2026.

9. KOPPULA, P.; ZHANG, Y.; ZHUANG, L.; GAN, B. Cystine transporter regulation of pentose phosphate pathway dependency and disulfide stress exposes a targetable metabolic vulnerability in cancer. **Nature Cell Biology**, v. 22, n. 4, p. 476-486, 2020. DOI: [10.1038/s41556-020-0496-x](https://doi.org/10.1038/s41556-020-0496-x). Disponível em: https://www.nature.com/articles/s41556-020-0496-x. Acesso em: 15 jun. 2026.

10. GOLDMAN, M. J.; CRAFT, B.; HASTIE, M.; REPEČKA, K.; MCDADE, F.; KAMATH, A.; BANERJEE, A.; LUO, Y.; ROGERS, D.; BROOKS, A. N.; ZHU, J.; HAUSSLER, D. Visualizing and interpreting cancer genomics data via the Xena platform. **Nature Biotechnology**, v. 38, n. 6, p. 675-678, 2020. DOI: [10.1038/s41587-020-0546-8](https://doi.org/10.1038/s41587-020-0546-8). Disponível em: https://www.nature.com/articles/s41587-020-0546-8. Acesso em: 20 jun. 2026.

11. GTEX CONSORTIUM. The Genotype-Tissue Expression (GTEx) project. **Nature Genetics**, v. 45, n. 6, p. 580-585, 2013. DOI: [10.1038/ng.2653](https://doi.org/10.1038/ng.2653). Disponível em: https://www.nature.com/articles/ng.2653. Acesso em: 28 jun. 2026.

12. RITCHIE, M. E.; PHIPSON, B.; WU, D.; HU, Y.; LAW, C. W.; SHI, W.; SMYTH, G. K. limma powers differential expression analyses for RNA-sequencing and microarray studies. **Nucleic Acids Research**, v. 43, n. 7, p. e47, 2015. DOI: [10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007). Disponível em: https://academic.oup.com/nar/article/43/7/e47/2414268. Acesso em: 4 jul. 2026.

13. VIVIAN, J.; RAO, A. A.; NOTHAFT, F. A.; KETCHUM, C.; ARMSTRONG, J.; NOVAK, A.; PFEIL, J.; NARKIZIAN, J.; DERAN, A. D.; MUSSELMAN-BROWN, A.; SCHMIDT, H.; AMSTUTZ, P.; CRAFT, B.; GOLDMAN, M.; ROSENBLOOM, K.; CLINE, M.; O'CONNOR, B.; HANNA, M.; BIRGER, C.; KENT, W. J.; PATTERSON, D. A.; JOSEPH, A. D.; ZHU, J.; ZARANEK, S.; GETZ, G.; HAUSSLER, D.; PATEN, B. Toil enables reproducible, open source, big biomedical data analyses. **Nature Biotechnology**, v. 35, n. 4, p. 314-316, 2017. DOI: [10.1038/nbt.3772](https://doi.org/10.1038/nbt.3772). Disponível em: https://www.nature.com/articles/nbt.3772. Acesso em: 8 jul. 2026.

14. WU, D.; SMYTH, G. K. Camera: a competitive gene set test accounting for inter-gene correlation. **Nucleic Acids Research**, v. 40, n. 17, p. e133, 2012. DOI: [10.1093/nar/gks461](https://doi.org/10.1093/nar/gks461). Disponível em: https://academic.oup.com/nar/article/40/17/e133/2411347. Acesso em: 11 jul. 2026.

15. UHLÉN, M.; FAGERBERG, L.; HALLSTRÖM, B. M.; LINDSKOG, C.; OKSVOLD, P.; MARDINOGLU, A.; SIVERTSSON, Å.; KAMPF, C.; SJÖSTEDT, E.; ASPLUND, A.; OLSSON, I.; EDLUND, K.; LUNDBERG, E.; NAVANI, S.; SZIGYARTO, C. A.; ODEBERG, J.; DJUREINOVIC, D.; TAKANEN, J. O.; HOBER, S.; ALM, T.; EDQVIST, P. H.; BERLING, H.; TEGEL, H.; MULDER, J.; ROCKBERG, J.; NILSSON, P.; SCHWENK, J. M.; HAMSTEN, M.; VON FEILITZEN, K.; FORSBERG, M.; PERSSON, L.; JOHANSSON, F.; ZWAHLEN, M.; VON HEIJNE, G.; NIELSEN, J.; PONTÉN, F. Tissue-based map of the human proteome. **Science**, v. 347, n. 6220, p. 1260419, 2015. DOI: [10.1126/science.1260419](https://doi.org/10.1126/science.1260419). Disponível em: https://www.science.org/doi/10.1126/science.1260419. Acesso em: 16 jul. 2026.

16. ARAN, D.; SIROTA, M.; BUTTE, A. J. Systematic pan-cancer analysis of tumour purity. **Nature Communications**, v. 6, p. 8971, 2015. DOI: [10.1038/ncomms9971](https://doi.org/10.1038/ncomms9971). Disponível em: https://www.nature.com/articles/ncomms9971. Acesso em: 20 jul. 2026.

---

## Licença

MIT License. Veja [LICENSE](LICENSE).

Copyright © 2026 Kamila da Conceição Loureiro, Ryan de Paulo Santos, Letícia Maria Dias Freitas, Ivine Souza Silva, Maria Eduarda Peixoto Soares Pecly.

---

*Última atualização: 12 de julho de 2026. Versão 3.1.0.*
