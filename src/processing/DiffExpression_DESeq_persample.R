library(DESeq2)
library(edgeR)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(UpSetR)

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
  str_replace("SLX-\\d{5} UDI\\d{4} MCF7 ", "") #%>%
  #str_extract("(.*?)(?=( ))")

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
  str_replace("SLX-\\d{5} UDI\\d{4} T47D ", "")# %>%
  #str_extract("(.*?)(?=( ))")

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

#Prepare to add int he DESeq2 design
coldata <- data.frame(row.names=colnames(final_mat), group, cell_line, shRNA)
#from coldata remove control rows 
col_data <- coldata[-grep("CTRL", rownames(coldata)), ]

#Add control row !!!!
tmp <- data.frame("CTRL", "CTRL", "CTRL")
rownames(tmp) <- "CTRL"
colnames(tmp) <- c("group", "cell_line", "shRNA")
col_data <- rbind(col_data, tmp)

# Grep CTRL columns
dat <- final_mat[,grep("CTRL", colnames(final_mat))]
head(dat)
# mean expression row-wise / u need integer cause deseq2 accepts raw counts 
#Which makes sense, because none of those values are integers! 
CTRL <- as.integer(rowMeans(dat))

# Drop CTRL columns
final_mat <- final_mat[, -grep("CTRL", colnames(final_mat))]
new_mat <- cbind(final_mat, as.data.frame(CTRL))

# <><><><><> Create DESeqDataSet <><><><><>

#If you want to use DESeq2 to analyze your data, you need to get count data which will be strictly integer values.
new_group <- group[-grep("CTRL", group)]
new_group <- c(new_group, "CTRL")

dds <- DESeqDataSetFromMatrix(countData=new_mat, colData=col_data, design = ~ group)
dds$group <- factor(dds$group, levels = unique(new_group))

dds <- DESeq(dds, parallel = TRUE)
#This object is saved in /output/TransformedData_DerivedData/DE_DESeq_bayes/deseq.Rdata


out_dir <- "/output/TransformedData_DerivedData/DE_DESeq_bayes/"
parallel::mclapply(unique(colData(dds)$group), function(g) {
  print(g)
  g <- as.character(g)
  tmp_res <-
    results(dds, pAdjustMethod = "BH", contrast = c("group", g, "CTRL"))
  resOrdered <- tmp_res[order(tmp_res$pvalue), ]
  write.csv(as.data.frame(resOrdered),
            file = paste0(out_dir, gsub(" ", "_", g), "_results_DESeq.csv")
  )
}, mc.cores = 18)

