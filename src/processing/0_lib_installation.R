#2020Sept15
#Anna Mathioudaki

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#RNAseq preprocessing & normalization
BiocManager::install("edgeR")
BiocManager::install("DESeq2") # we dont use this 

install.packages("ggplot2")

# Useful visualization libraries ~ not necessary but useful
install.packages("pheatmap")
BiocManager::install("EnhancedVolcano")
install.packages("UpSetR")
