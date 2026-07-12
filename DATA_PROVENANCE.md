# Data Provenance

## Primary Data Source

| Propriedade | Valor |
|-------------|-------|
| **Plataforma** | [UCSC Xena](https://xenabrowser.net/) |
| **Dataset** | `TcgaTargetGtex_RSEM_Hugo_norm_count` |
| **Hub** | `toil-xena-hub` (pipeline Toil RNA-seq recompute) |
| **URL de download** | `https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz` |
| **Data de acesso** | 11 de julho de 2026 |
| **Formato** | TSV (compressed GZIP ~1.2 GB; descompressed ~6 GB) |

## Arquivos de dados

### kidney.tsv (versionado — ~220 KB)

Arquivo reduzido contendo apenas os genes das vias metabólicas (66 colunas de genes) para todas as 445 amostras renais. Usado como referência para identificar as amostras renais durante o download do transcriptoma completo.

| Propriedade | Valor |
|-------------|-------|
| SHA256 | `29154696504ce365e681d9c319fe352a6c84c5ae87798c8f1cce3c11159f7ea2` |
| Amostras | 445 |
| Colunas de genes | 66 (genes do metabolismo central) |

### kidney_transcriptome.tsv (NÃO versionado — ~178 MB)

Transcriptoma completo: todos os 58.581 genes para as 445 amostras renais. Gerado por `scripts/download_full_transcriptome.py`.

| Propriedade | Valor |
|-------------|-------|
| Amostras | 445 |
| Genes (brutos) | 58.581 |
| Genes (após filtragem) | 31.633 |
| Genes removidos (baixa expressão) | 26.948 |
| Genes removidos (símbolos ambíguos) | consultar pipeline_v3.R |

### Outros arquivos em data/raw/

| Arquivo | Descrição | Versionado? |
|---------|-----------|:-----------:|
| `kidney.tsv.gz` | kidney.tsv comprimido | Não |
| `kidney_full.tsv` | Transcriptoma sem filtragem | Não |
| `kidney_full.tsv.gz` | Comprimido | Não |
| `kidney_tcga_kirp.tsv.gz` | Subset TCGA-KIRP apenas | Não |

## Composição das Amostras

| Grupo | n | Study | Sample Type | Código |
|-------|---|-------|-------------|--------|
| KIRP Tumor | 288 | TCGA | Primary Tumor | `-01` |
| GTEx Normal | 28 | GTEx | Normal Tissue (Kidney Cortex) | prefixo `GTEX-` |
| TCGA Normal (KIRP-adjacente) | 32 | TCGA | Solid Tissue Normal | `-11`, mesmo participante que tumor KIRP |
| TCGA Normal (outros projetos) | 97 | TCGA | Solid Tissue Normal | `-11`, projetos KIRC, KICH, outros |

**Total: 445 amostras**

## Processamento dos Dados

### Pipeline Toil (pré-processamento pelo UCSC Xena)

1. Dados brutos de RNA-seq (HTSeq counts) do TCGA e GTEx
2. Processados pelo pipeline Toil RNA-seq recompute (VIVIAN et al., 2017)
3. Normalização: RSEM (RNA-Seq by Expectation-Maximization)
4. Transformação: log₂(norm_count + 1) aplicada pelo UCSC Xena
5. Integração: TCGA + GTEx combinados em matriz única

### Filtragem local (pipeline_v3.R)

1. **Símbolos ambíguos:** genes com `?` no símbolo removidos
2. **Baixa expressão:** genes com expressão > 1 em ≤ 10% das 445 amostras removidos
3. **Resultado:** 58.581 → 31.633 genes

### Escala de expressão

| Propriedade | Valor |
|-------------|-------|
| **Escala** | log₂(norm_count + 1) |
| **Mínimo** | 0 |
| **Máximo** | 20,1 |
| **Mediana** | 11,15 |
| **% zeros** | 9,6% |
| **Valores ausentes** | 0 |

> ⚠️ **Importante:** A matriz já é fornecida em escala log₂ pelo UCSC Xena. NENHUMA transformação logarítmica adicional foi aplicada neste estudo. Aplicar log₂ novamente produziria valores incorretos.

## Metadados das Amostras

| Coluna | Descrição |
|--------|-----------|
| `sample` | ID da amostra no UCSC Xena |
| `study` | TCGA ou GTEx |
| `sample_type` | Primary_Tumor, Solid_Tissue_Normal, Normal_Tissue |
| `condition` | KIRP, TCGA_Normal, Normal_GTEx |
| `participant` | ID do participante (3 primeiros segmentos do código TCGA) |
| `paired_normal` | TRUE se for normal adjacente de participante KIRP |

## Genes das Vias Metabólicas

| Via | KEGG ID | Genes no KEGG | Genes na Matriz |
|-----|---------|:-------------:|:---------------:|
| Glicólise / Gliconeogênese | hsa00010 | 67 | 64 |
| Via das Pentoses Fosfato | hsa00030 | 31 | 30 |
| Ciclo do Ácido Cítrico (TCA) | hsa00020 | 30 | 29 |
| **União (genes únicos)** | — | **110** | **106** |

**Genes ausentes (n = 4):** *G6PC1*, *PRPS1L1*, *RPEL1*, *SUCLA2* — removidos por baixa expressão.

Mapeamento congelado em: `results/tables/pathway_gene_membership.csv` (KEGG release 119.0, julho/2026).

## Caveats Importantes

1. **Confundimento condição-coorte (GTEx):** Todos os tumores KIRP são do TCGA; todos os normais GTEx são de autópsias. O efeito da condição (tumor vs. normal) e o efeito da coorte (TCGA vs. GTEx) são perfeitamente confundidos. A análise KIRP vs. GTEx é estritamente exploratória.

2. **Tecido adjacente não é necessariamente normal:** O tecido adjacente ao tumor pode apresentar alterações de campo, inflamação e contaminação tumoral.

3. **Normais TCGA de projetos mistos:** Os 97 normais TCGA não-KIRP provêm de múltiplos projetos (KIRC, KICH, outros contextos renais). Não representam um grupo homogêneo.

4. **Gene-level, não isoform-level:** Os dados medem expressão gênica total. Inferências sobre isoformas específicas não são suportadas.

5. **Dados pré-transformados:** A escala log₂(norm_count + 1) é fornecida pelo UCSC Xena. Interpretações sobre contagens absolutas requerem transformação reversa (2^x − 1), que é aproximada.

## Referência do Dataset

GOLDMAN, M. J. et al. Visualizing and interpreting cancer genomics data via the Xena platform. **Nature Biotechnology**, v. 38, n. 6, p. 675-678, 2020. DOI: [10.1038/s41587-020-0546-8](https://doi.org/10.1038/s41587-020-0546-8).

VIVIAN, J. et al. Toil enables reproducible, open source, big biomedical data analyses. **Nature Biotechnology**, v. 35, n. 4, p. 314-316, 2017. DOI: [10.1038/nbt.3772](https://doi.org/10.1038/nbt.3772).
