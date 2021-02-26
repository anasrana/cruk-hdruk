library(DESeq2)
library(ggplot2)
library(tidyverse)

# The here command allows you to pass simpler paths based on project root in R
mcf_file <- "/data/turing_data_dump/MCF7/Gene/count_matrix.csv"
t47d_file <- "/data/turing_data_dump/T47D/Gene/count_matrix.csv"

# Read the data fromt eh cell lines and extract some information from col names
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
for (i in colnames(mcf7_mat)){
  gene <- str_replace(i, "SLX-\\d{5} UDI\\d{4} MCF7 ", "") %>% str_extract("(.*?)(?=( ))")
  tmp_sh <- str_replace(i, paste0("SLX-\\d{5} UDI\\d{4} MCF7 ", gene), "") %>% str_extract("[^ ]+")
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
for (i in colnames(t47d_mat)){
  gene <- str_replace(i, "SLX-\\d{5} UDI\\d{4} T47D ", "") %>% str_extract("(.*?)(?=( ))")
  tmp_sh <- str_replace(i, paste0("SLX-\\d{5} UDI\\d{4} T47D ", gene), "") %>% str_extract("[^ ]+")
  shRNA2 <- c(shRNA2, tmp_sh)
}

#Merge matrices coming from different cell lines 
final_mat <- cbind(mcf7_mat, t47d_mat)
group <- c(group1, group2)
cell_line <- c(cell_line1, cell_line2)
shRNA <- c(shRNA1, shRNA2)

#mcf7_mat
coldata1 <- data.frame(row.names=colnames(mcf7_mat), group1, cell_line1, shRNA1)
dds1 <- DESeqDataSetFromMatrix(countData=mcf7_mat, colData=coldata1, design = ~ group1)
dds1$group1 <- factor(dds1$group1, levels = unique(group1))
dds1 <- DESeq(dds1)
rld1 <- rlogTransformation(dds1, blind=TRUE)

pdf("PCA_rlog_KO_mcf7.pdf", height = 5, width = 5)
plotPCA(rld1, intgroup= "group1") + geom_point(size = 0) + theme(text=element_text(size=15)) + ggtitle("MCF7 cell line")
dev.off()

#Save dds1 and rld1 objects
saveRDS(dds1, file = "/output/TransformedData_DerivedData/DESeq_Robjects/DESeq_obj_mcf7.rds")
saveRDS(rld1, file = "/output/TransformedData_DerivedData/DESeq_Robjects/rlog_counts_mcf7.rds")

#t47d_mat
coldata2 <- data.frame(row.names=colnames(t47d_mat), group2, cell_line2, shRNA2)
dds2 <- DESeqDataSetFromMatrix(countData=t47d_mat, colData=coldata2, design = ~ group2)
dds2$group2 <- factor(dds2$group2, levels = unique(group2))
dds2 <- DESeq(dds2)
rld2 <- rlogTransformation(dds2, blind=TRUE)

pdf("PCA_rlog_KO_t47d.pdf", height = 5, width = 5)
plotPCA(rld2, intgroup= "group2") + geom_point(size = 0) + theme(text=element_text(size=15)) + ggtitle("T47D cell line")
dev.off()

#Save dds1 and rld1 objects
saveRDS(dds2, file = "/output/TransformedData_DerivedData/DESeq_Robjects/DESeq_obj_t47d.rds")
saveRDS(rld2, file = "/output/TransformedData_DerivedData/DESeq_Robjects/rlog_counts_t47d.rds")
