# PATHWAY EXPANSION RATIONALE — Central Carbon Metabolism in KIRP

**Date:** 2026-07-11  
**Version:** v2.0.0 preregistration  
**Status:** Committed before observing pathway results

---

## 1. Rationale for hsa00010 — Glycolysis / Gluconeogenesis

The glycolytic pathway (KEGG hsa00010) is the central energy-yielding pathway of carbohydrate metabolism and was the focus of v1. Numerous studies document altered glycolytic gene expression in renal cell carcinoma (RCC), driven by HIF-1α stabilization in clear cell RCC [PMID:34071836] and by distinct mechanisms in papillary RCC, including MET mutations and fumarate accumulation [PMID:38253811, PMID:35216667]. The Warburg effect — maintenance of high glycolytic rate even under aerobic conditions — is a hallmark of cancer metabolism [PMID:36330953]. However, the transcriptomic profile of glycolysis/gluconeogenesis genes in KIRP specifically remains underexplored relative to KIRC [PMID:23792563].

**Key references:**
- Cancer Cell Metabolism in Hypoxia (PMID:34071836, OA-PMC, Int J Mol Sci 2021)
- Metabolic alterations in hereditary/sporadic RCC (PMID:38253811, OA-PMC, Nat Rev Nephrol 2024)
- Fumarate inhibits PTEN in type 2 papillary RCC (PMID:35216667, Mol Cell 2022)
- Tumor glycolysis review (PMID:36330953, Semin Cancer Biol 2022)

---

## 2. Rationale for hsa00030 — Pentose Phosphate Pathway

The pentose phosphate pathway (PPP, KEGG hsa00030) branches from glycolysis at glucose-6-phosphate and serves two critical functions: (i) production of NADPH for reductive biosynthesis and antioxidant defense, and (ii) production of ribose-5-phosphate for nucleotide biosynthesis [PMID:25037503]. Cancer cells upregulate PPP to support proliferation and manage oxidative stress [PMID:37612403].

In RCC, the PPP has been implicated in redox homeostasis. G6PD, the rate-limiting enzyme of PPP, is relevant to kidney biology [PMID:27755120]. The cystine/glutamate antiporter xCT regulates PPP dependency in cancer [PMID:32231310]. However, the transcriptomic status of PPP genes in KIRP remains largely uncharacterized. Given that KIRP type 2 is associated with fumarate hydratase (FH) deficiency and consequent oxidative stress [PMID:35216667], the PPP may be particularly relevant for NADPH production in this subtype.

**Key references:**
- The pentose phosphate pathway and cancer (PMID:25037503, OA-PMC, Trends Biochem Sci 2014)
- The pentose phosphate pathway in health and disease (PMID:37612403, OA-PMC, Nat Metab 2023)
- G6PD and the kidney (PMID:27755120, Curr Opin Nephrol Hypertens 2017)
- Cystine transporter regulation of PPP (PMID:32231310, OA-PMC, Nat Cell Biol 2020)

---

## 3. Rationale for hsa00020 — Citrate Cycle (TCA Cycle)

The tricarboxylic acid (TCA) cycle (KEGG hsa00020) occupies the central hub of oxidative metabolism, coupling carbohydrate, lipid, and amino acid catabolism to ATP production via oxidative phosphorylation, while providing biosynthetic intermediates (citrate for fatty acid synthesis, α-ketoglutarate for amino acid synthesis) [PMID:31032474].

In RCC, the TCA cycle is profoundly affected by subtype-specific mutations. KIRC is defined by VHL loss leading to pseudo-hypoxia, while a subset of KIRP (type 2) harbors FH mutations causing fumarate accumulation and succination of proteins [PMID:38253811, PMID:35216667]. Complex I of the mitochondrial electron transport chain has been recently implicated in kidney cancer metastasis [PMID:39143213]. Targeting TCA-deficient RCC is an active therapeutic area [PMID:36128328].

The interplay between glycolysis, PPP, and TCA cycle — collectively the "central carbon metabolism" — is a systems-level property [PMID:28502706]. Disruption in one pathway propagates to others through shared intermediates: glucose-6-phosphate partitions between glycolysis and PPP; pyruvate feeds the TCA cycle or is fermented to lactate; citrate exported from mitochondria supports lipid synthesis.

**Key references:**
- Krebs cycle metabolites in immunity and cancer (PMID:31032474, OA-PMC, Nat Metab 2019)
- Mitochondrial complex I promotes kidney cancer metastasis (PMID:39143213, OA-PMC, Nature 2024)
- Targeting Krebs-cycle-deficient RCC (PMID:36128328, OA-PMC, Oncotarget 2022)
- Reprogramming of central carbon metabolism in cancer stem cells (PMID:28502706, BBA Mol Basis Dis 2017)

---

## 4. Metabolic Interconnection of the Three Pathways

```
Glucose
  |
  v
Glucose-6-P ----(G6PD)----> 6-Phosphogluconolactone ----> PPP (hsa00030)
  |                                                         |
  | (GPI)                                              NADPH + R5P
  v
Fructose-6-P
  |
  | (PFK)
  v
GLYCOLYSIS (hsa00010) ------> Pyruvate
                                    |
                          (PDH)     |     (PC)
                               v    |    v
                              Acetyl-CoA   Oxaloacetate
                                   |          |
                                   v          v
                              TCA CYCLE (hsa00020)
                                   |
                              NADH, FADH2, GTP
                                   |
                              Citrate → Fatty acid synthesis
                              α-KG → Amino acid synthesis
                              Succinate → HIF stabilization
                              Fumarate → Succination (KIRP type 2)
```

The three pathways are not independent modules but an integrated metabolic network. Gene expression changes in one pathway may be compensated or exacerbated by changes in another.

---

## 5. Relevance to KIRP

KIRP is genetically and metabolically distinct from KIRC. While KIRC is driven by VHL/HIF/pseudo-hypoxia, KIRP subtypes include:
- **Type 1:** MET alterations, gains of chromosomes 7 and 17
- **Type 2:** CDKN2A silencing, SETD2 mutations, FH deficiency in a subset

FH-deficient KIRP type 2 accumulates fumarate, a TCA cycle intermediate, which stabilizes HIF-1α through competitive inhibition of α-KG-dependent dioxygenases and promotes succination of cellular proteins including PTEN [PMID:35216667]. This directly connects the TCA cycle to tumor biology in papillary RCC.

The pentose phosphate pathway may be relevant for KIRP through its role in NADPH production for managing oxidative stress, which may be elevated in FH-deficient tumors.

---

## 6. Knowledge Gaps

1. The transcriptomic profile of PPP genes in KIRP has not been systematically characterized alongside glycolysis
2. The integrated analysis of all three central carbon metabolism pathways in KIRP has not been reported
3. The robustness of transcriptomic findings across different comparator strategies (TCGA-normal, GTEx-normal) has not been evaluated for these three pathways
4. The implications of these integrative findings for precision health have not been articulated

---

## 7. Testable Hypotheses (preregistered)

1. KIRP tumors will show transcriptomic alterations in genes of glycolysis/gluconeogenesis, PPP, and TCA cycle compared to normal kidney tissue
2. The degree of alteration will differ across the three pathways, reflecting their distinct metabolic roles
3. Some genes will show robust alterations across comparator strategies; others will be comparator-sensitive
4. The interconnection nodes (G6PD, PKM, PDHA1, PDHB, PCK1, PCK2) will show coordinated alterations

---

## 8. Limitations (preregistered)

1. RNA expression does not measure metabolic flux or enzyme activity
2. Bulk RNA-seq does not resolve cell-type heterogeneity
3. The TCGA KIRP dataset does not routinely distinguish type 1 from type 2 papillary RCC
4. Pathway membership is based on KEGG annotations which may not capture all relevant genes
5. Transcript-level (isoform) resolution is not available for PKM and other alternatively spliced genes

---

*Document committed before observing any v2 pathway results.*
