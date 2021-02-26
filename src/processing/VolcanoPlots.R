library(EnhancedVolcano)

list_files <- list.files("/output/TransformedData_DerivedData/DE_DESeq/")
out_dir <- "/output/Images/InitialAnalysis/VolcanoPlots/"

for (f in list_files){
  if (f != "README"){
    gene <- strsplit(f, split = "_")[[1]][3]
    DE_all <- read.csv(paste0("/scratch/shared/DE_DESeq/", f), row.names = 1)
    DE_all$gene <- rownames(DE_all)
    DE_all$baseMean <- NULL
    DE_all$lfcSE <- NULL
    DE_all$stat <- NULL
    colnames(DE_all) <- c("log2FoldChange", "p_val", "pvalue", "gene")
    
    de_all <- as.data.frame(DE_all)
    
    res2 <- de_all
    keyvals <- ifelse(
      res2$log2FoldChange < -1, '#F8766D',
      ifelse(res2$log2FoldChange > 1, '#00BFC4',
             'grey30'))
    keyvals[is.na(keyvals)] <- 'black'
    names(keyvals)[keyvals == '#F8766D'] <- 'control'
    names(keyvals)[keyvals == '#00BFC4'] <- 'shRNA'
    names(keyvals)[keyvals == 'grey30'] <- 'p-value only'
    
    pdf(paste0(out_dir, gene, "_VolcanoPlot.pdf"), width = 6, height = 5)
    p <- EnhancedVolcano(res2, 
                         lab = res2$gene, 
                         x = "log2FoldChange", 
                         y = "pvalue", 
                         labSize = 3.0, 
                         FCcutoff = 1, 
                         ylab = "-log10(pval_adj)", 
                         xlab = "log2FC", 
                         title = gene, 
                         col = c("grey30", "grey30", "grey30", "#F8766D"), 
                         subtitle = "", 
                         colAlpha = 1/2, 
                         colCustom = keyvals, 
                         axisLabSize = 15, 
                         titleLabSize = 15, 
                         subtitleLabSize = 10, 
                         caption = NULL)
    plot(p)
    dev.off()
  }
}
