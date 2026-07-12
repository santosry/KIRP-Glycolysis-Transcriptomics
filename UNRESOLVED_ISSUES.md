# UNRESOLVED_ISSUES.md

## Pendências que requerem reexecução do pipeline

### 1. Transcriptoma completo vs painel direcionado
- **Arquivo:** `data/raw/kidney_transcriptome.tsv` (186 MB, 58.581 genes)
- **Ação:** Reexecutar `01_load_data.R` usando esta matriz, filtrar para ~31.633 genes, rodar limma no universo completo, extrair os 106 genes a posteriori
- **Impacto:** Permite camera, GSEA e ORA com universos estatisticamente válidos

### 2. Exclusão de normais TCGA não-KIRP da filtragem
- **Ação:** Aplicar filtro de baixa expressão apenas sobre amostras KIRP (288) + KIRP-adjacentes (32) + GTEx (28) = 348 amostras
- **Impacto:** Evita remoção de genes com expressão específica em KIRP

### 3. Rede de coexpressão sobre amostras mistas
- **Ação:** Calcular correlações usando expressão residualizada pelo modelo pareado, ou separadamente por condição
- **Impacto:** Remove correlações espúrias entre tumores e normais

### 4. Ambiente de execução
- **Problema:** R segfault no Git Bash (Msys2) ao carregar certos pacotes
- **Workaround:** Usar PowerShell para scripts R
- **Impacto:** Renderização do PDF requer PowerShell; `run_all.R` funciona via PowerShell
