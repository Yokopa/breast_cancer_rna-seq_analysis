# 467713-HS2023: RNA-sequencing - Master in Bioinformatics and Computational Biology, University of Bern
# breast_cancer_rna-seq_analysis  

# Load libraries
library(DESeq2)
library(ggplot2)
library(dplyr)

# Read the data into a data frame
# -------------------------------
counts_data <- read.table("featurecounts.Rmatrix.txt", header = TRUE, sep = "\t", quote = "")
head(counts_data)

# Create a DataFrame with metadata
# --------------------------------
sample_names <- colnames(counts_data)[-1] # Extract sample names from column headers

colData <- data.frame(
  sampleName = sample_names,
  group = c(rep("HER2", 3), rep("NonTNBC", 3), rep("Normal", 3), rep("TNBC", 3))
)

# Construct DESEQDataSet Object
# -----------------------------
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~ group,
                              tidy = TRUE)

# Run DESEQ function
# ------------------
dds <- DESeq(dds)
res <- results(dds)

# Look at the results table
summary(res)

# Plot PCA
# --------
# Perform variance stabilizing transformation (VST)
vsdata <- vst(dds, blind=TRUE)
plotPCA(vsdata, intgroup="group")

res <- as.data.frame((res))
total_genes <- sum(res$padj < 0.05, na.rm=TRUE)
# Count upregulated and downregulated genes
upregulated_genes <- sum(res$log2FoldChange > 0 & res$padj < 0.05, na.rm = TRUE)
downregulated_genes <- sum(res$log2FoldChange < 0 & res$padj < 0.05, na.rm = TRUE)

# Print the counts
cat("Tital genes:", total_genes, "\n")
cat("Up-regulated genes:", upregulated_genes, "\n")
cat("Down-regulated genes:", downregulated_genes, "\n")

# Create data frames of 3 selected contrasts
# ------------------------------------------
# TNBC vs HER2
TNBC_HER2 <- results(dds, contrast = c("group", "TNBC", "HER2"))
TNBC_HER2 <- as.data.frame(TNBC_HER2)
# TNBC vs NonTNBC
TNBC_NonTNBC <- results(dds, contrast = c("group", "TNBC", "NonTNBC"))
TNBC_NonTNBC <- as.data.frame(TNBC_NonTNBC)
# TNBC vs Normal
TNBC_Normal <- results(dds, contrast = c("group", "TNBC", "Normal"))
TNBC_Normal <- as.data.frame(TNBC_Normal)
# Her2 vs Normal
#HER2_Normal <- results(dds, contrast = c("group", "HER2", "Normal"))
#HER2_Normal <- as.data.frame(TNBC_Normal)
# NonTNBC vs Normal
#NonTNBC_Normal <- results(dds, contrast = c("group", "NonTNBC", "Normal"))
#NonTNBC_Normal <- as.data.frame(TNBC_Normal)

# Save the results
#write.csv(TNBC_HER2, "deseq2_TNBC_HER2.all.csv")
#write.csv(TNBC_NonTNBC, "deseq2_TNBC_NonTNBC.all.csv")
#write.csv(TNBC_Normal, "deseq2_TNBC_Normal.all.csv")


total_genes <- sum(TNBC_NonTNBC$padj < 0.05, na.rm=TRUE)
# Count upregulated and downregulated genes
upregulated_genes <- sum(TNBC_NonTNBC$log2FoldChange > 0 & TNBC_NonTNBC$padj < 0.05, na.rm = TRUE)
downregulated_genes <- sum(TNBC_NonTNBC$log2FoldChange < 0 & TNBC_NonTNBC$padj < 0.05, na.rm = TRUE)

# Print the counts
cat("TNBC vs NonTNBC - Total genes:", total_genes, "\n")
cat("TNBC vs NonTNBC - Up-regulated genes:", upregulated_genes, "\n")
cat("TNBC vs NonTNBC - Down-regulated genes:", downregulated_genes, "\n")


# Results filtering
# -----------------
# Sort DataFrame by adjusted p-value
TNBC_HER2 <- TNBC_HER2[order(TNBC_HER2$padj),] 
TNBC_NonTNBC <- TNBC_NonTNBC[order(TNBC_NonTNBC$padj),]
TNBC_Normal <- TNBC_Normal[order(TNBC_Normal$padj),]
#HER2_Normal <- HER2_Normal[order(HER2_Normal$padj),]
#NonTNBC_Normal <- NonTNBC_Normal[order(NonTNBC_Normal$padj),] 

# Filter based on p adjusted value
filtered_TNBC_HER2 <- TNBC_HER2 %>% filter(TNBC_HER2$padj < 0.05)
filtered_TNBC_NonTNBC <- TNBC_NonTNBC %>% filter(TNBC_NonTNBC$padj < 0.05)
filtered_TNBC_Normal <- TNBC_Normal %>% filter(TNBC_Normal$padj < 0.05)
#filtered_HER2_Normal <- HER2_Normal %>% filter(HER2_Normal$padj < 0.05)
#filtered_NonTNBC_Normal <- NonTNBC_Normal %>% filter(NonTNBC_Normal$padj < 0.05)

# Look at the results
head(filtered_TNBC_HER2)
head(filtered_TNBC_NonTNBC)
head(filtered_TNBC_Normal)

# Translate gene ID's to gene names
# ---------------------------------
library(org.Hs.eg.db)

# TNBC vs HER2
filtered_TNBC_HER2$symbol <- mapIds(org.Hs.eg.db, keys = rownames(filtered_TNBC_HER2), keytype = "ENSEMBL", column = "SYMBOL")
#TNBC vs NonTNBC
filtered_TNBC_NonTNBC$symbol <- mapIds(org.Hs.eg.db, keys = rownames(filtered_TNBC_NonTNBC), keytype = "ENSEMBL", column = "SYMBOL")
# TNBC vs Normal
filtered_TNBC_Normal$symbol <- mapIds(org.Hs.eg.db, keys = rownames(filtered_TNBC_Normal), keytype = "ENSEMBL", column = "SYMBOL")
# TNBC vs HER2
#filtered_HER2_Normal$symbol <- mapIds(org.Hs.eg.db, keys = rownames(filtered_HER2_Normal), keytype = "ENSEMBL", column = "SYMBOL")
# NonTNBC vs Normal
filtered_NonTNBC_Normal$symbol <- mapIds(org.Hs.eg.db, keys = rownames(filtered_NonTNBC_Normal), keytype = "ENSEMBL", column = "SYMBOL")

# Look at the results
head(filtered_TNBC_HER2)
head(filtered_TNBC_NonTNBC)
head(filtered_TNBC_Normal)

# Save the results
#write.csv(filtered_TNBC_HER2, "deseq2_TNBC_HER2.filtered.csv")
write.csv(filtered_TNBC_NonTNBC, "deseq2_TNBC_NonTNBC.filtered.csv")
#write.csv(filtered_TNBC_Normal, "deseq2_TNBC_Normal.filtered.csv")

# Results visualization with Volcano Plot
# ---------------------------------------
library(EnhancedVolcano)
library(ggrepel)

# TNBC vs HER2
EnhancedVolcano(filtered_TNBC_HER2,
                lab = filtered_TNBC_HER2$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'TNBC vs HER2',
                pointSize = 1.0,
                labSize = 4.0,
                colAlpha = 0.4,
                legendPosition = 'None',
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                max.overlaps = 20)

# TNBC vs NonTNBC
EnhancedVolcano(filtered_TNBC_NonTNBC,
                lab = filtered_TNBC_NonTNBC$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'TNBC vs NonTNBC',
                pointSize = 1.0,
                labSize = 4.0,
                colAlpha = 0.4,
                legendPosition = 'None',
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                max.overlaps = 20)

# TNBC vs Normal
EnhancedVolcano(filtered_TNBC_Normal,
                lab = filtered_TNBC_Normal$symbol,
                x = 'log2FoldChange',
                y = 'padj',
                title = 'TNBC vs Normal',
                pointSize = 1.0,
                labSize = 3.0,
                colAlpha = 0.4,
                legendPosition = 'None',
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                max.overlaps = 20)

# TNBC vs NonTNBC with customized colors for the report "
keyvals <- ifelse(
  filtered_TNBC_NonTNBC$log2FoldChange < 0 & filtered_TNBC_NonTNBC$padj < 0.05, 'steelblue3',
  ifelse(filtered_TNBC_NonTNBC$log2FoldChange > 0 & filtered_TNBC_NonTNBC$padj < 0.05, 'palevioletred', 'grey'))

names(keyvals)[keyvals == 'steelblue3'] <- 'Down-regulated'
names(keyvals)[keyvals == 'palevioletred'] <- 'Up-regulated'
names(keyvals)[keyvals == 'grey'] <- 'Not significant'

EnhancedVolcano(
  filtered_TNBC_NonTNBC,
  lab = filtered_TNBC_NonTNBC$symbol,
  x = 'log2FoldChange',
  y = 'padj',
  title = 'TNBC vs NonTNBC',
  pCutoff = 0.05, FCcutoff = 1,
  #pointSize = 3.0,
  #labSize = 6.0,
  colAlpha = 0.3,
  #drawConnectors = TRUE,
  #widthConnectors = 0.2,
  max.overlaps = 20,
  colCustom = keyvals)

# Save plot
ggsave("TNBC_NonTNBC_volcanoplot.png",
       plot = last_plot(),
       scale = 1,
       )

# Investigation of 3 selected genes: ROPN1B, ROPN and ACTG2
# ---------------------------------------------------------
# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Select genes of interest (in order: ELF5, ROPN1B, ROPN and ACTG2)
genes_of_interest <- c("ENSG00000135374", "ENSG00000114547", "ENSG00000065371", "ENSG00000163017") 

# Check if the selected genes are present in the dataset
selected_genes_counts <- normalized_counts[rownames(normalized_counts) %in% genes_of_interest, ]

# Print the selected genes and check if data is available: ENSG00000135374 = ELF5, ENSG00000114547 = ROPN1b, ENSG00000065371 = ROPN, ENSG00000163017 = ACTG2
print(selected_genes_counts)

# Over-representation analysis with clusterProfiler
# -------------------------------------------------
library(clusterProfiler)
library(AnnotationDbi)

# Get the gene set for DE genes
de_genes_TNBC_HER2 <- rownames(filtered_TNBC_HER2)
de_genes_TNBC_NonTNBC <- rownames(filtered_TNBC_NonTNBC)
de_genes_TNBC_Normal <- rownames(filtered_TNBC_Normal)

# Get the universe set of genes
universe_genes <- counts_data$Geneid

# Run enrichGO for TNBC vs HER2
enrichment_TNBC_HER2 <- enrichGO(gene = de_genes_TNBC_HER2,
                                 universe = universe_genes,
                                 OrgDb = org.Hs.eg.db,
                                 keyType = "ENSEMBL",
                                 ont = "BP",
                                 pvalueCutoff = 0.05,
                                 readable = TRUE)
# Run enrichGO for TNBC vs NonTNBC
enrichment_TNBC_NonTNBC <- enrichGO(gene = de_genes_TNBC_NonTNBC,
                                 universe = universe_genes,
                                 OrgDb = org.Hs.eg.db,
                                 keyType = "ENSEMBL",
                                 ont = "BP",
                                 pvalueCutoff = 0.05,
                                 readable = TRUE)
# Run enrichGO for TNBC vs Normal
enrichment_TNBC_Normal <- enrichGO(gene = de_genes_TNBC_Normal,
                                    universe = universe_genes,
                                    OrgDb = org.Hs.eg.db,
                                    keyType = "ENSEMBL",
                                    ont = "BP",
                                    pvalueCutoff = 0.05,
                                    readable = TRUE)

# View the results
as.data.frame(enrichment_TNBC_HER2)
as.data.frame(enrichment_TNBC_NonTNBC)
as.data.frame(enrichment_TNBC_Normal)

# Visualization of functional enrichment result
# ---------------------------------------------
#library(enrichplot)

# Barplots

# TNBC vs HER2
bar_enr.TNBC_HER2 <- plot(barplot(enrichment_TNBC_HER2, showCategory = 10, title = "TNBC vs HER2"))

# Save barplot
ggsave("enrichment_TNBC_HER2_barplot.png",
       plot = last_plot(),
       scale = 1,
)

# TNBC vs NonTNBC
bar_enr.TNBC_NonTNBC <- plot(barplot(enrichment_TNBC_NonTNBC, showCategory = 30, title = "TNBC vs NonTNBC"))
bar_enr.TNBC_NonTNBC <- plot(barplot(enrichment_TNBC_NonTNBC, showCategory = 10))

# Save barplot
ggsave("enrichment_TNBC_NonTNBC_barplot.png",
       plot = last_plot(),
       scale = 1,
)

# TNBC vs Normal
bar_enr.TNBC_Normal <- plot(barplot(enrichment_TNBC_Normal, showCategory = 20, title = "TNBC vs Normal"))

# Save barplot
ggsave("enrichment_TNBC_Normal_barplot.png",
       plot = last_plot(),
       scale = 1,
)

# Dotplot

# TNBC vs HER2
plot(dotplot(enrichment_TNBC_HER2))

# TNBC vs NonTNBC
plot(dotplot(enrichment_TNBC_NonTNBC, showCategory = 15))

# Save dotplot
ggsave("enrichment_TNBC_NonTNBC_dotplot.png",
       plot = last_plot(),
       scale = 1,
)

# TNBC vs Normal
plot(dotplot(enrichment_TNBC_Normal))

# Cnetplot 

# TNBC vs HER2
cnetplot(enrichment_TNBC_HER2, circular = TRUE, colorEdge = TRUE)
# TNBC vs NonTNBC
cnetplot(enrichment_TNBC_NonTNBC, circular = TRUE, colorEdge = TRUE)