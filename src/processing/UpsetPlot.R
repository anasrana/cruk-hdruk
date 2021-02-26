library(UpSetR)

#List of silenced genes
group <- c("SUMO3", "NCOA3", "SUMO1", "TFAP2C", "CREBBP", "NRIP1", "GRHL2", "NR2F2", "ZMIZ1", "TRIM33", "SUMO2", "RARA", "TRIM24", "FOXA1", "EP300", "DPF2", "GATA3", "TRIM28", "ESRRA")

#Directory which contains the differential expression DESeq2 output from /DSGSept2020CRUKCI/src/processing/DiffExpression_DESeq.R
out_dir <- "/output/TransformedData_DerivedData/DE_DESeq/"

# Summary of the DE genes across different knock outs
sig_genes <- c()
for (gene in unique(group)){
    tmp <- read.csv(paste0(out_dir, "CTRL_vs_", gene, "_results_DESeq.csv"), row.names = 1)
    tmp <- na.omit(tmp)
    tmp_sig <- rownames(tmp[tmp$padj < 0.05, ])
    sig_genes[[gene]] <- tmp_sig
    print(paste0(gene, " has ", length(tmp_sig), " significantly differentially expressed"))
}

pdf("/output/Images/InitialAnalysis/Upset_DEgenes.pdf", width = 7, height = 7)
upset(fromList(sig_genes),
      nsets = 20,
      matrix.color = "darkred")
dev.off()
