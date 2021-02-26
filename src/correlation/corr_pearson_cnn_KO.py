# Script to compute pairwise pearson correlations to be used in CNN
# Knockout genes are added if not present in the selection list

import pandas as pd
from scipy.stats import pearsonr
import numpy as np
import pickle as pkl
from multiprocessing import Pool
from itertools import combinations

def calc_pearson(data, gene, gene_list):
    pearson = pd.DataFrame(0, columns=gene_list, index=gene_list)
    for i, j in combinations(gene_list, 2):
        pearson.loc[i,j] = pearsonr(data.loc[i], data.loc[j])[0]

    pearson += pearson.T
    np.fill_diagonal(pearson.values, 1.0)
    pkl.dump(pearson, open(f"../../TransformedData_DerivedData/correlation_network/corr_coeff_pearson_cnn_{gene}_KO.pkl", "wb"))

if __name__ == "__main__":
    df = pd.read_csv("../../TransformedData_DerivedData/norm_count.csv", index_col = 0)
    # Multiindex columns to be (cellline, gene, shRNA, replicate)
    df.columns=pd.MultiIndex.from_tuples([tuple(x.split()[-4:]) for x in df.columns])
    systems = sorted(set(x[1] for x in df.columns))
    genes = sorted(set(x[1] for x in df.columns))
    genes.remove("CTRL")
    # CNN gene list supplied by Jessica Dafflon
    with open("../../TransformedData_DerivedData/CNN/histogram/gene_list_cnn.csv") as f:
        large_genes = f.read().strip().split("\n")
    # Throw in the knockout genes
    large_genes += genes
    large_genes = sorted(set(large_genes))

    with Pool(20) as p:
        p.starmap(calc_pearson, [(df.xs(g, axis=1, level=1, drop_level=False).loc[large_genes], g, large_genes) for g in systems])