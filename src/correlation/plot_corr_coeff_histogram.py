# Plot the histograms for correlation coefficients

import pandas as pd
import pickle as pkl
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

if __name__ == "__main__":
    systems = pkl.load(open("../../TransformedData_DerivedData/correlation_network/corr_coeff_20.pkl", "rb")).keys()
    # Pearson
    data = []
    for i in systems:
        corr = pkl.load(open(f"../../TransformedData_DerivedData/correlation_network/corr_coeff_pearson_{i}_KO.pkl", "rb"))
        data.append(corr.values.flatten())
    plt.figure(figsize=(10, 5))
    sns.histplot(np.concatenate(data))
    plt.xlabel("Pearson correlation coefficient")
    plt.savefig("../../Images/correlation_network/pearson_histogram.pdf")

    # Spearman
    data = []
    for i in systems:
        corr = pkl.load(open(f"../../TransformedData_DerivedData/correlation_network/corr_coeff_spearman_{i}_KO.pkl", "rb"))
        data.append(corr.values.flatten())
    plt.figure(figsize=(10, 5))
    sns.histplot(np.concatenate(data))
    plt.xlabel("Pearson correlation coefficient")
    plt.savefig("../../Images/correlation_network/spearman_histogram.pdf")
    