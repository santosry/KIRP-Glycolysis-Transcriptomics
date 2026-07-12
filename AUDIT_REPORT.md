# AUDIT_REPORT.md — Auditoria Científica v3.2.0

**Data:** 12 de julho de 2026  
**Repositório:** `santosry/KIRP-Glycolysis-Transcriptomics`  
**Commit:** `3c15477`

---

## 🔴 CRÍTICO — Universo gênico: contradição fundamental

**Evidência:** `data/raw/kidney_transcriptome.tsv` (186 MB, 58.581 genes × 445 amostras) existe e está versionado. Contudo, o pipeline atual (`01_load_data.R`) utiliza `data/raw/kidney.tsv` (336 KB, 109 genes).

**Impacto:** Todas as análises de enriquecimento (camera, GSEA, ORA) operam sobre um universo de ~106 genes metabólicos, tornando seus resultados parcialmente circulares. O manuscrito descreve o estudo como análise "transcriptômica" quando na prática é um painel gênico direcionado.

**Correção necessária:**
1. O pipeline deve ser reexecutado sobre `kidney_transcriptome.tsv` com limma nos 31.633 genes pós-filtragem, extraindo os 106 genes a posteriori.
2. OU o manuscrito deve ser reclassificado como "análise de painel gênico direcionado", removendo todas as alegações de análise transcriptômica global.

**Status:** 🔴 NÃO RESOLVIDO — requer reexecução do pipeline completo.

---

## 🔴 CRÍTICO — "Transcriptoma" no título

**Evidência:** O título diz "Transcriptoma das vias glicolítica...". Uma análise de 106 genes não constitui análise transcriptômica.

**Correção:** Alterar para "Perfil de expressão gênica das vias glicolítica, das pentoses fosfato e do ciclo do ácido cítrico no carcinoma renal papilar: análise pareada e sensibilidade ao tecido de referência".

**Status:** ✅ Corrigido no manuscrito.

---

## 🔴 CRÍTICO — Camera com universo inválido

**Evidência:** `camera()` foi executado sobre 106 genes. O teste competitivo requer milhares de genes de background.

**Correção:** Remover camera das conclusões. Mantê-lo apenas como registro de transparência com forte ressalva.

**Status:** ✅ Corrigido no manuscrito.

---

## 🔴 CRÍTICO — GSEA com universo inválido

**Evidência:** `GSEA_KEGG.csv` contém 352 gene sets testados sobre ranking de 106 genes. O `p.adjust` para hsa01200 = 0.012 e hsa00010 = 0.031 são reportados como significativos, mas o universo é insuficiente.

**Correção:** Remover GSEA das conclusões. Mencionar apenas como análise exploratória.

**Status:** ✅ Corrigido no manuscrito.

---

## 🟡 ALTO — Amostras TCGA de outros projetos na filtragem

**Evidência:** 129 normais TCGA incluem 32 KIRP-adjacentes + 97 de outros projetos (KIRC, KICH). A filtragem de baixa expressão (`rowMeans > 1 > 0.1`) foi aplicada sobre todas as 445 amostras, incluindo os 97 normais não-KIRP.

**Impacto:** Genes com expressão específica em KIRP podem ter sido removidos se ausentes nos normais de outros projetos. A filtragem deveria usar apenas amostras KIRP e KIRP-adjacentes.

**Correção:** Refazer filtragem usando apenas amostras KIRP (n=288) + KIRP-adjacentes (n=32) + GTEx (n=28) = 348 amostras, excluindo os 97 normais de outros projetos.

**Status:** 🟡 PENDENTE — requer reexecução do pipeline.

---

## 🟡 ALTO — "Extração a posteriori" inapropriada

**Evidência:** O manuscrito afirma que os genes foram "extraídos a posteriori do modelo transcriptômico único". Na prática, a matriz foi obtida por seleção direcionada no UCSC Xena.

**Correção:** Remover a expressão. Descrever honestamente como "painel gênico pré-especificado das três vias KEGG".

**Status:** ✅ Corrigido no manuscrito.

---

## 🟡 ALTO — ORA circular

**Evidência:** A ORA usa os 106 genes como universo. Enriquecimento de "glicólise" nesse universo é tautológico.

**Correção:** Adicionar ressalva explícita de circularidade. Remover interpretações causais.

**Status:** ✅ Corrigido no manuscrito.

---

## 🟡 MÉDIO — Rede de coexpressão sobre amostras mistas

**Evidência:** `04_ppi_network.R` calcula correlação de Pearson sobre todas as 445 amostras (tumores + normais misturados), gerando correlações artificiais pela separação entre condições.

**Correção:** Usar expressão residualizada ou calcular redes separadas.

**Status:** 🟡 PENDENTE — requer modificação do script.

---

## 🟡 MÉDIO — Número de genes inconsistente

**Evidência:** O pipeline reporta 106 genes finais (109 raw - 3 removidos). O manuscrito menciona 106. O README menciona 106. A S1 tem 106 genes. Consistente entre si, mas não reflete que G6PC1 está ausente do Xena e 18 genes apareciam duplicados (deduplicados manualmente).

**Status:** ✅ Números reconciliados. Processo de deduplicação documentado.

---

## 🟢 BAIXO — Referências com datas de acesso futuras

**Evidência:** Datas de acesso variam de 15 mar. 2026 a 20 jul. 2026. A data real de acesso ao Xena foi 11-12 jul. 2026.

**Correção:** Uniformizar para datas reais de acesso.

**Status:** ✅ Ajustado para 11-12 jul. 2026.

---

## 🟢 BAIXO — Figuras com numeração inconsistente

**Evidência:** O manuscrito refere-se a "Tabela 3" para a tabela de concordância, mas há múltiplas tabelas com numeração sobreposta.

**Status:** ✅ Corrigido.

---

## RESUMO

| Gravidade | Total | Resolvidos | Pendentes |
|-----------|:-----:|:----------:|:---------:|
| 🔴 Crítico | 4 | 3 | 1 |
| 🟡 Alto | 3 | 1 | 2 |
| 🟡 Médio | 2 | 1 | 1 |
| 🟢 Baixo | 2 | 2 | 0 |

**Pendências que requerem reexecução do pipeline:**
1. Usar `kidney_transcriptome.tsv` (58K genes) como matriz base
2. Excluir 97 normais TCGA não-KIRP da filtragem
3. Corrigir rede de coexpressão (residualizar ou separar por condição)
