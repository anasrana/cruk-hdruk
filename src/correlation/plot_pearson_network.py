import pandas as pd
import networkx as nx
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt
from multiprocessing import Pool

def normalise(p, n):
    # normalise positve and negative correlation coefficients
    p = np.array(p)
    n = np.array(n)
    small = min(min(p), min(abs(n)))
    large = max(max(p), max(abs(n)))
    p = (p - small) / large
    n = (n + small) / large
    return p, n

def corr_to_g(knockout, cut=0, norm=False):
    # convert correlation dataframe to a graph
    df = pkl.load(open(f"../../TransformedData_DerivedData/correlation_network/corr_coeff_pearson_{knockout}_KO.pkl", "rb"))
    G = nx.Graph()
    pcolors = []
    ncolors = []
    positive = []
    negative = []
    genes = [x for x in systems if x != "CTRL"]
    node_colours = np.array(abs(df[abs(df) > cut])[genes].sum())
    for x, y in np.argwhere(abs(np.triu(df.values)) > cut):
        v1 = df.columns[x]
        v2 = df.columns[y]
        if v1 in genes and v2 in genes:
            if df.loc[v1, v2] > 0:
                G.add_edge(v1, v2, weight=df.loc[v1, v2])
                pcolors.append(df.loc[v1, v2])
                positive.append((v1, v2))
            elif df.loc[v1, v2] < 0:
                G.add_edge(v1, v2, weight=abs(df.loc[v1, v2]))
                ncolors.append(abs(df.loc[v1, v2]))
                negative.append((v1, v2))
    fig = plt.figure(figsize=(10, 10))
    pos = nx.circular_layout(sorted(genes))
    if norm:
        pcolors, ncolors = normalise(pcolors, ncolors)
    nx.draw_networkx_nodes(G, pos, alpha=1.0, node_size=2000, node_color=node_colours, cmap=plt.cm.Greens, nodelist=genes)
    nx.draw_networkx_edges(G, pos, edge_color=pcolors, edgelist= positive, edge_cmap=plt.cm.Reds, edge_vmin=0.0, width=2.0, alpha=0.4, min_source_margin=25, min_target_margin=25)
    nx.draw_networkx_edges(G, pos, edge_color=ncolors, edgelist= negative, edge_cmap=plt.cm.Blues, edge_vmin=0.0, width=2.0, alpha=0.4, min_source_margin=25, min_target_margin=25)
    nx.draw_networkx_labels(G, pos, font_size=10)
    if knockout == "CTRL":
        plt.title("Pearson correlation network of the control sample", fontsize=20)
    else:
        plt.title(f"Pearson correlation network of the {knockout} knockout experiment", fontsize=20)
    plt.tight_layout()
    fig.savefig(f"../../Images/correlation_network/correlation_network_pearson_{knockout}.pdf")

if __name__ == "__main__":
    # get the names of systems from the correlation file
    systems = pkl.load(open("../../TransformedData_DerivedData/correlation_network/corr_coeff_20.pkl", "rb")).keys()
    # alternatively...
    # systems = ['CREBBP','CTRL','DPF2','EP300','ESRRA','FOXA1','GATA3','GRHL2','NCOA3',
    #             'NR2F2','NRIP1','RARA','SUMO1','SUMO2','SUMO3','TFAP2C','TRIM24','TRIM28','TRIM33','ZMIZ1']
    with Pool(20) as p:
        p.starmap(corr_to_g, [(s, 0.8, True) for s in systems])