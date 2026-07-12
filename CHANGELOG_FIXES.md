# CHANGELOG_FIXES.md

## v3.3.0 — Auditoria Científica (12 jul 2026)

### Correções aplicadas ao manuscrito

1. **Título:** "Transcriptoma" → "Perfil de expressão gênica" (análise restrita a 106 genes, não transcriptoma completo)
2. **Removido "extração a posteriori":** O estudo usa painel gênico pré-especificado, não extração do transcriptoma completo
3. **Camera:** Removido das conclusões; mantido com ressalva de universo inválido
4. **GSEA:** Reclassificado como exploratório; removido de conclusões
5. **ORA:** Adicionada ressalva de circularidade (universo de genes metabólicos)
6. **Terminologia:** "Up/Down" → "expressão aumentada/reduzida"; "tecido normal" → "tecido não tumoral adjacente"
7. **7 genes compartilhados:** Corrigido "duas ou três vias" → "duas vias" (todos pertencem a exatamente 2)
8. **Número de genes:** Uniformizado para 106 em todo o texto
9. **Referências:** Datas de acesso ajustadas para 11-12 jul 2026
10. **Abstract removido** (não necessário para esta versão)
11. **Agradecimentos removidos** (solicitado pelo autor)

### Pendências não resolvidas (requerem reexecução do pipeline)

- Usar `kidney_transcriptome.tsv` como matriz base (58K genes)
- Excluir 97 normais TCGA não-KIRP da filtragem
- Corrigir rede de coexpressão (residualização por condição)
