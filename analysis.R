# =============================================================
# Cervical Cancer RNA-seq Analysis
# Dataset: GSE168009
# Comparison: NDB (resistant) vs DCB (sensitive)
# =============================================================

library(GEOquery)
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)

# =============================================================
# STEP 1: Load Data
# =============================================================

# Load metadata
gse <- getGEO("GSE168009", GSEMatrix = TRUE)

# Extract phenotype data (sample information)
pheno <- pData(phenoData(gse[[1]]))

# Load raw counts
counts_file <- "GSE168009/GSE168009_Raw_count.txt.gz"
counts <- read.delim(counts_file, row.names = 1)

# Quick look at the data
dim(counts)        # how many genes x samples?
head(counts)       # first few rows
colnames(counts)   # sample names

# =============================================================
# STEP 2: Create metadata table
# =============================================================

# Create a dataframe describing each sample
col_data <- data.frame(
  sample = colnames(counts),
  condition = c("NDB", "NDB", "NDB", "NDB", "DCB", "DCB", "DCB", "DCB", "DCB"),
  row.names = colnames(counts)
)

col_data$condition <- factor(col_data$condition, levels = c("DCB", "NDB"))

# Check it looks right
col_data

# =============================================================
# STEP 3: Create DESeq2 object
# =============================================================

# Build the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = col_data,
  design = ~ condition
)

# Filter out lowly expressed genes
# Keep only genes with at least 10 reads total across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Check how many genes remain after filtering
dim(dds)

# =============================================================
# STEP 4: Run DESeq2
# =============================================================

dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("condition", "NDB", "DCB"))

# Summary of results
summary(res)

# Look at the top results
head(res[order(res$padj), ])

# =============================================================
# STEP 5: Visualisation — Volcano Plot
# =============================================================

# Convert results to dataframe for plotting
res_df <- as.data.frame(res)

# Volcano plot
EnhancedVolcano(res_df,
                lab = rownames(res_df),
                x = "log2FoldChange",
                y = "padj",
                title = "NDB vs DCB — Cervical Cancer Chemoradiotherapy Resistance",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0
)

# =============================================================
# STEP 6: Pathway Analysis
# =============================================================

# Install required packages
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library(clusterProfiler)
library(org.Hs.eg.db)

# Get significant genes (padj < 0.05, |log2FC| > 1)
sig_genes <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)

# Convert gene symbols to Entrez IDs
gene_ids <- bitr(rownames(sig_genes),
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

# Run GO enrichment analysis
go_results <- enrichGO(gene = gene_ids$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       readable = TRUE)

# Visualise
dotplot(go_results, showCategory = 15)

# =============================================================
# STEP 7: Save outputs
# =============================================================

# Create output directory
dir.create("outputs", showWarnings = FALSE)

# Save DESeq2 results as CSV
write.csv(as.data.frame(res[order(res$padj), ]),
          file = "outputs/DESeq2_results.csv")

# Save significant genes only
write.csv(sig_genes[order(sig_genes$padj), ],
          file = "outputs/significant_genes.csv")

# Save volcano plot
png("outputs/volcano_plot.png", width = 10, height = 8, units = "in", res = 300)
EnhancedVolcano(res_df,
                lab = rownames(res_df),
                x = "log2FoldChange",
                y = "padj",
                title = "NDB vs DCB — Cervical Cancer Chemoradiotherapy Resistance",
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0)
dev.off()

# Save dotplot
png("outputs/GO_dotplot.png", width = 10, height = 8, units = "in", res = 300)
dotplot(go_results, showCategory = 15)
dev.off()

# =============================================================
# STEP 8: PCA Plot
# =============================================================

# Variance stabilising transformation — needed for PCA
vsd <- vst(dds, blind = FALSE)

# Plot PCA
plotPCA(vsd, intgroup = "condition") +
  theme_minimal() +
  ggtitle("PCA — NDB vs DCB") +
  scale_color_manual(values = c("DCB" = "#2166ac", "NDB" = "#d6604d"))

# =============================================================
# STEP 9: Heatmap
# =============================================================

BiocManager::install("pheatmap")
library(pheatmap)

# Get top 30 most significant genes
top_genes <- head(order(res$padj), 30)

# Extract their normalised counts
mat <- assay(vsd)[top_genes, ]

# Scale by row (so we're comparing relative expression, not absolute)
mat_scaled <- t(scale(t(mat)))

# Annotation for columns (samples)
annotation <- data.frame(
  condition = col_data$condition,
  row.names = colnames(mat)
)

# Plot
pheatmap(mat_scaled,
         annotation_col = annotation,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("#2166ac", "white", "#d6604d"))(100),
         main = "Top 30 DE Genes — NDB vs DCB")
png("outputs/heatmap_top30.png", width = 10, height = 12, units = "in", res = 300)
pheatmap(mat_scaled,
         annotation_col = annotation,
         show_rownames = TRUE,
         show_colnames = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("#2166ac", "white", "#d6604d"))(100),
         main = "Top 30 DE Genes — NDB vs DCB")
dev.off()
