# loading library for simple manipulations and reading data
library(tidyverse)

# The here command allows you to pass simpler paths based on project root in R
mcf_file <- here::here("data/MCF7/Gene/count_matrix.csv")
t47d_file <- here::here("data/T47D/Gene/count_matrix.csv")

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

mcf_df %>% distinct(k_out)

t47d_df <- read_csv(t47d_file) %>%
    pivot_longer(-X1, names_to = "sample", values_to = "count") %>%
    rename(genes = X1) %>%
    mutate(
        s_id = str_extract(sample, "UDI\\d{4}"),
        sample = str_replace(sample, "SLX-\\d{5} UDI\\d{4} T47D ", ""),
        k_out = str_extract(sample, "(.*?)(?=( ))"),
        c_line = "T47D"
    )

# Create simple boxplot of raw data
b_plot <- bind_rows(mcf_df, t47d_df) %>%
    ggplot(aes(x = sample, y = count)) +
    geom_boxplot() +
    facet_grid(c_line ~ k_out, scales = "free_x", space = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave("<WHERE TO SAVE>", b_plot, width = 20, height = 10)

# ==============================================================================
# Prepare data for edgeR analysis
# ==============================================================================

library(edgeR)
library(RColorBrewer)
library(mixOmics)


tmp_df <-
    read_csv(mcf_file) %>%
    dplyr::filter(!str_detect(X1, "NA\\.[0-9]{0,3}"))

mcf7_mat <-
    tmp_df[, -1] %>%
    as.matrix() %>%
    magrittr::set_rownames(tmp_df$X1)

group <-
    colnames(mcf7_mat) %>%
    str_replace("SLX-\\d{5} UDI\\d{4} MCF7 ", "") %>%
    str_extract("(.*?)(?=( ))")

# Create a list object used for edgeR
mcf7_dat <- DGEList(counts = mcf7_mat, group = group)

# calculate normalisation factor and dispersion of genes
mcf7_dat <- calcNormFactors(mcf7_dat, method = "TMM")

mcf7_dat <- estimateGLMRobustDisp(mcf7_dat, verbose = T)

# Info on normalisation can be found looking at
head(mcf7_dat$samples)

# extract normalised counts
norm_count <- log2(cpm(mcf7_dat, normalized.lib.sizes = TRUE, log = FALSE) + 1)

# Histogram of first sample in norm count
hist(norm_count[, 1])

# http://www.nathalievialaneix.eu/doc/html/TP1_normalization.html

plotMDS(mcf7_dat,
    top = 1000, labels = mcf7_dat$samples$group,
    col = as.numeric(mcf7_dat$samples$group), cex = 2
)


# http://www.nathalievialaneix.eu/doc/html/TP1_normalization.html

# create sample distance matrix
sample_dists <- as.matrix(dist(t(norm_count)))

# create heatmap of distance matrix
cim_col <- colorRampPalette(rev(brewer.pal(9, "Reds")))(16)
cim(sample_dists,
    color = cim_col,
    symkey = FALSE, row.cex = 0.7, col.cex = 0.7
)
