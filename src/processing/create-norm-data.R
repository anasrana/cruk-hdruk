# loading library for simple manipulations and reading data
library(tidyverse)
library(edgeR)

# The here command allows you to pass simpler paths based on project root in R
mcf_file <- here::here("/data/turing_data_dump/MCF7/Gene/count_matrix.csv")
t47d_file <- here::here("/data/turing_data_dump/T47D/Gene/count_matrix.csv")

tmp_df <-
    full_join(
        read_csv(mcf_file),
        read_csv(t47d_file),
        by = "X1"
    ) %>%
    dplyr::filter(!str_detect(X1, "NA\\.[0-9]{0,3}"))

count_mat <-
    tmp_df[, -1] %>%
    as.matrix() %>%
    magrittr::set_rownames(tmp_df$X1)

group <-
    colnames(count_mat) %>%
    str_replace("SLX-\\d{5} UDI\\d{4} ", "") %>%
    str_extract("(MCF7|T47D) (.*?)(?=( ))")

# Create a list object used for edgeR
all_dat <- DGEList(counts = count_mat, group = group, remove.zeros = T)

# calculate normalisation factor and dispersion of genes
all_dat <- calcNormFactors(all_dat, method = "TMM")

all_dat <- estimateGLMRobustDisp(all_dat, verbose = T)

# extract normalised counts
norm_count_l <-
    log2(cpm(all_dat, normalized.lib.sizes = TRUE, log = FALSE) + 1) %>%
    as_tibble(rownames = "gene")

norm_count <- cpm(all_dat, normalized.lib.sizes = TRUE, log = FALSE) %>%
    as_tibble(rownames = "gene")


write_csv(norm_count, "norm_count.csv")
write_csv(norm_count_l, "norm_count-log.csv")
