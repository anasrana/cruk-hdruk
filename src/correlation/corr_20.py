import pandas as pd
from itertools import combinations
from scipy.stats import pearsonr
import numpy as np
import pickle as pkl

if __name__ == "__main__":
    df = pd.read_csv("../../TransformedData_DerivedData/norm_count.csv", index_col = 0)
    # Multiindex columns to be (cellline, gene, shRNA, replicate)
    df.columns=pd.MultiIndex.from_tuples([tuple(x.split()[-4:]) for x in df.columns])

    # get all knockout genes
    genes = sorted(set(x[1] for x in df.columns))
    genes.remove("CTRL")

    # all perturbations + control
    systems = sorted(set(x[1] for x in df.columns))

    corr = {x: None for x in systems}
    for s in systems:
        # cross section of a system and keep rows of genes only
        data = df.xs(s, axis=1, level=1, drop_level=False).loc[genes]
        # initialise empty dataframe
        pearson = pd.DataFrame(0, columns=genes, index=genes)
        # fill upper triangle
        for i, j in combinations(genes, 2):
            pearson.loc[i,j] = pearsonr(data.loc[i], data.loc[j])[0]
        # copy over to lower triangle and fill diagonal
        pearson += pearson.T
        np.fill_diagonal(pearson.values, 1.0)
    
        corr[s] = pearson
    pkl.dump(corr, open("../../TransformedData_DerivedData/correlation_network/corr_coeff_20.pkl", "wb"))