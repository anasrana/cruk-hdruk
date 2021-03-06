import pandas as pd
from scipy.stats import kendalltau
import numpy as np
import pickle as pkl
from multiprocessing import Pool
from itertools import combinations

def calc_kendall(data, gene, gene_list):
    kendall = pd.DataFrame(0, columns=gene_list, index=gene_list)
    for i, j in combinations(gene_list, 2):
        try:
            kendall.loc[i,j] = kendalltau(data.loc[i], data.loc[j])[0]
        except:
            print(i, j, gene)
            exit()

    kendall += kendall.T
    np.fill_diagonal(kendall.values, 1.0)
    pkl.dump(kendall, open(f"aboli_kendall_corr_coeff_large_{gene}.pkl", "wb"))

if __name__ == "__main__":
    df = pd.read_csv("norm_count.csv", index_col = 0)
    # Multiindex columns to be (cellline, gene, shRNA, replicate)
    df.columns=pd.MultiIndex.from_tuples([tuple(x.split()[-4:]) for x in df.columns])
    systems = sorted(set(x[1] for x in df.columns))
    corr = {x: None for x in systems}
    genes = sorted(set(x[1] for x in df.columns))
    genes.remove("CTRL")
    large_genes = []
    for g in genes:
        data = pd.read_csv(f"/shared/DE_DESeq/CTRL_vs_{g}_results_DESeq.csv", index_col=0)
        large_genes += list(data[data.padj < 0.05].sort_values("padj").index[:100])
    large_genes = sorted(set(large_genes))
    with Pool(10) as p:
        p.starmap(calc_kendall, [(df.xs(g, axis=1, level=1, drop_level=False).loc[large_genes], g, large_genes) for g in systems])