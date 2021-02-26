# 2020Sept16
# AnnaMathioudaki

# Load Libraries
pacman::p_load(DESeq2, tidyverse, pheatmap)

# Define output directory
out_dir <- "/output/TransformedData_DerivedData/DE_DESeq/"

# The here command allows you to pass simpler paths based on project root in R
# mcf_file <- "/data/turing_data_dump/MCF7/Gene/count_matrix.csv"
# t47d_file <- "/data/turing_data_dump/T47D/Gene/count_matrix.csv"
mcf_file <- "data/MCF7/Gene/count_matrix.csv"
###
# Read the data from the cell lines and extract some information from col names
mcf_df <- read_csv(mcf_file) %>%
    pivot_longer(-X1, names_to = "sample", values_to = "count") %>%
    rename(genes = X1) %>%
    mutate(
        s_id = str_extract(sample, "UDI\\d{4}"),
        sample = str_replace(sample, "SLX-\\d{5} UDI\\d{4} MCF7 ", ""),
        k_out = str_extract(sample, "(.*?)(?=( ))"),
        c_line = "MCF7"
    )

t47d_df <- read_csv(t47d_file) %>%
    pivot_longer(-X1, names_to = "sample", values_to = "count") %>%
    rename(genes = X1) %>%
    mutate(
        s_id = str_extract(sample, "UDI\\d{4}"),
        sample = str_replace(sample, "SLX-\\d{5} UDI\\d{4} T47D ", ""),
        k_out = str_extract(sample, "(.*?)(?=( ))"),
        c_line = "T47D"
    )

###
tmp_df <-
    read_csv(mcf_file) %>%
    dplyr::filter(!str_detect(X1, "NA\\.[0-9]{0,3}"))

mcf7_mat <-
    tmp_df[, -1] %>%
    as.matrix() %>%
    magrittr::set_rownames(tmp_df$X1)

group1 <-
    colnames(mcf7_mat) %>%
    str_replace("SLX-\\d{5} UDI\\d{4} MCF7 ", "") %>%
    str_extract("(.*?)(?=( ))")

cell_line1 <- rep("mcf7", length(group1))

shRNA1 <- c()
for (i in colnames(mcf7_mat)) {
    gene <-
        str_replace(i, "SLX-\\d{5} UDI\\d{4} MCF7 ", "") %>%
        str_extract("(.*?)(?=( ))")
    tmp_sh <-
        str_replace(i, paste0("SLX-\\d{5} UDI\\d{4} MCF7 ", gene), "") %>%
        str_extract("[^ ]+")

    shRNA1 <- c(shRNA1, tmp_sh)
}


tmp_df <-
    read_csv(t47d_file) %>%
    dplyr::filter(!str_detect(X1, "NA\\.[0-9]{0,3}"))

t47d_mat <-
    tmp_df[, -1] %>%
    as.matrix() %>%
    magrittr::set_rownames(tmp_df$X1)

group2 <-
    colnames(t47d_mat) %>%
    str_replace("SLX-\\d{5} UDI\\d{4} T47D ", "") %>%
    str_extract("(.*?)(?=( ))")

cell_line2 <- rep("t47d", length(group2))

shRNA2 <- c()
for (i in colnames(t47d_mat)) {
    gene <-
        str_replace(i, "SLX-\\d{5} UDI\\d{4} T47D ", "") %>%
        str_extract("(.*?)(?=( ))")
    tmp_sh <-
        str_replace(i, paste0("SLX-\\d{5} UDI\\d{4} T47D ", gene), "") %>%
        str_extract("[^ ]+")

    shRNA2 <- c(shRNA2, tmp_sh)
}


# Merge matrices coming from different cell lines
final_mat <- cbind(mcf7_mat, t47d_mat)
group <- c(group1, group2)
cell_line <- c(cell_line1, cell_line2)
shRNA <- c(shRNA1, shRNA2)

# Prepare to add info in the DESeq2 design
coldata <- data.frame(row.names = colnames(final_mat), group, cell_line, shRNA)
head(coldata)

# <><><><><> Create DESeqDataSet <><><><><>
dds <-
    DESeqDataSetFromMatrix(
        countData = final_mat, colData = coldata,
        design = ~ group + cell_line
    )
dds$group <- factor(dds$group, levels = unique(group))
dds <- DESeq(dds)

# And Differential expression
for (gene in unique(group)) {
    if (gene != "CTRL") {
        tmp_res <-
            results(
                dds,
                pAdjustMethod = "BH",
                contrast = c("group", gene, "CTRL")
            )
        resOrdered <- tmp_res[order(tmp_res$pvalue), ]
        write.csv(
            as.data.frame(resOrdered),
            file = paste0(out_dir, "CTRL_vs_", gene, "_results_DESeq.csv")
        )
    }
}

# Rlog transform the data
rld <- rlogTransformation(dds, blind = TRUE)

pdf("/output/Images/InitialAnalysis/PCA_rlog/PCA_rlog_cellline.pdf",
    height = 5.3, width = 5
)

plotPCA(rld, intgroup = "cell_line") +
    geom_point(size = 0) +
    theme(text = element_text(size = 15)) +
    ggtitle("Both cell lines")

dev.off()

# Heatmap of log2FC for silenced genes across all perurbations
df2 <- as.data.frame(matrix(nrow = 20))

for (gene in unique(group)) {
    if (gene != "CTRL") {
        tmp <-
            read.csv(
                paste0(out_dir, "CTRL_vs_", gene, "_results_DESeq.csv"),
                row.names = 1
            )
        tmp_df <-
            data.frame(
                row.names = unique(group),
                tmp[unique(group), "log2FoldChange"]
            )
        colnames(tmp_df) <- gene
        df2 <- cbind(df2, tmp_df)
    }
}
df2$V1 <- NULL

pheatmap(
    df2,
    show_colnames = TRUE,
    show_rownames = TRUE, cluster_cols = FALSE, cluster_rows = FALSE
)

# Save dds and rld objects
saveRDS(dds,
    file = "/output/TransformedData_DerivedData/DESeq_Robjects/DESeq_obj_mcf7_t47d.rds"
)
saveRDS(rld,
    file = "/output/TransformedData_DerivedData/DESeq_Robjects/rlog_counts_mcf7_t47d.rds"
)