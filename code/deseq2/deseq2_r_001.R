#Installation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("pheatmap")
BiocManager::install("vsn")

#Input data from htseq
directory <- "/Users/emblapersson/Downloads/htseq/"

sampleFiles <- grep("htseq",list.files(directory),value=TRUE)
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = c("hp126", "hp126", "hp126", "r7", "r7", "r7"))
sampleTable$condition <- factor(sampleTable$condition)

#Load libraries
library("DESeq2")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")

#Differential expression analysis
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq

dds <- DESeq(ddsHTSeq)
res <- results(dds)
res <- res[!(res$pvalue>0.05 | is.na(res$pvalue)),]
res <- res[complete.cases(res), ]
res

significant_gene_names <- rownames(res)  
dds_filtered <- dds[significant_gene_names, ]

#p-values
resOrdered <- res[order(res$pvalue),]

#Log-fold change
resultsNames(dds_filtered)
resLFC <- lfcShrink(dds_filtered, coef="condition_r7_vs_hp126", type="apeglm")
resLFC

#Exploring and exporting results

#Volcano plot:
res$significant <- (abs(res$log2FoldChange) > 1 & -log10(res$pvalue) > 1.3)

# Create the volcano plot with similar criteria as in create_volcano_plot
volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), alpha = 0.5, size = 1.5) +  # Match alpha and size from create_volcano_plot
  scale_color_manual(values = c("gray", "red"), labels = c("Not Significant", "Significant")) +  # Color coding
  theme_minimal() +  # Minimalist plot style
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +  # Horizontal line at p = 0.05
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +  # Vertical lines for fold change thresholds
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-log10(p-value)",
    color = "Significance"
  )

volcano_plot

#MA-plot
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

# Heatmap
ntd <- normTransform(dds_filtered)
log2FC_values <- res$log2FoldChange
select <- order(abs(log2FC_values),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_filtered))
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

# Heatmap genes 21-40
select <- order(abs(log2FC_values),
                decreasing=TRUE)[21:40]
df <- as.data.frame(colData(dds_filtered))
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

# PCA plot
vst <- vst(dds)
plotPCA(vst)

vst_filtered <- vst(dds_filtered)
plotPCA(vst_filtered)
