# Cervical Cancer RNA-seq Analysis
## Transcriptomic signatures of chemoradiotherapy resistance in cervical cancer

---

## Overview
Differential gene expression analysis of a publicly available cervical cancer 
dataset comparing patients with no durable benefit (NDB, resistant) vs durable 
clinical benefit (DCB, sensitive) following platinum-based chemoradiotherapy.

**Dataset:** GSE168009 (GEO, NCBI)  
**Samples:** 9 patient tissue biopsies (4 NDB, 5 DCB)  
**Platform:** Illumina HiSeq 2500  

---

## Biological Question
Which genes and pathways are differentially expressed in cervical cancer patients 
who fail to respond to platinum-based chemoradiotherapy?

---

## Workflow
1. Data retrieval from GEO using GEOquery
2. Quality filtering (genes with ≥10 total reads retained)
3. Differential expression analysis with DESeq2
4. Visualisation — volcano plot
5. Pathway enrichment analysis — GO Biological Process (clusterProfiler)

---

## Key Findings
- 345 differentially expressed genes identified (padj < 0.05, |log2FC| > 1)
- Top upregulated in resistant: SLCO1B1, PRAME
- Top downregulated in resistant: OTOP3, CNTNAP3P2, DYNAP
- Pathway analysis reveals enrichment of extracellular matrix organisation 
  and immune signalling pathways in resistant tumours

---

## Tools & Packages
- R 4.6.0
- DESeq2
- clusterProfiler
- EnhancedVolcano
- GEOquery
- tidyverse

---

## Outputs
- `outputs/DESeq2_results.csv` — full results table
- `outputs/significant_genes.csv` — filtered significant genes
- `outputs/volcano_plot.png` — differential expression visualisation
- `outputs/GO_dotplot.png` — pathway enrichment visualisation

