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
