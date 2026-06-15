# Methods Summary

Este projeto analisa genes da via glicólise/gliconeogênese em carcinoma renal papilar (KIRP) usando dados transcriptômicos públicos TCGA/GTEx.

1. A via `hsa00010` foi obtida no KEGG com `KEGGREST`.
2. Símbolos gênicos humanos foram extraídos do campo `GENE`.
3. A matriz `kidney.tsv` foi importada de `data/raw/`.
4. Amostras foram classificadas como `Normal` ou `KIRP`.
5. A matriz de expressão foi convertida para formato numérico e transposta para genes em linhas e amostras em colunas.
6. A expressão diferencial foi estimada por `limma`, usando contraste `KIRP_vs_Normal`.
7. A correção de múltiplos testes usou Benjamini-Hochberg.
8. DEGs foram definidos por FDR < 0,05 e `|logFC| > 1`.
9. O volcano plot foi salvo em PNG com 300 dpi.
10. A rede PPI foi construída com STRING, filtrada por `combined_score >= 700` e salva em PNG com 300 dpi.

Análises de sobrevida, Kaplan-Meier, enriquecimento funcional, UpSet, dotplot e PPI 3D não fazem parte do fluxo final.
