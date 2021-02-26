import pandas as pd
import networkx as nx
import numpy as np
import pickle as pkl
from matplotlib import pyplot as plt
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
    G = nx.Graph()
    # inner = list(abs(df).sum().sort_values(ascending=False)[:10].index)
    df = pkl.load(open(f"../../TransformedData_DerivedData/correlation_network/corr_coeff_spearman_{knockout}_KO.pkl", "rb"))
    centrality = pkl.load(open(f"../../TransformedData_DerivedData/correlation_network/betweenness_centrality_spearman_{knockout}.pkl", "rb"))
    inner = sorted(centrality, key = lambda x: centrality[x], reverse=True)[:10]
    systems = pkl.load(open("../../TransformedData_DerivedData/correlation_network/corr_coeff_20.pkl", "rb")).keys()
    outter = [x for x in systems if x != "CTRL"]
    inner = [x for x in inner if x not in outter]
    shells = [inner, sorted(outter)]
    pcolors = []
    ncolors = []
    positive = []
    negative = []
    node_colours = []
    genes = inner + outter
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
    for c in G.nodes:
        if c in inner:
            node_colours.append("indianred")
        elif c in outter:
            node_colours.append("g")
    fig = plt.figure(figsize=(15, 15))
    pos = nx.shell_layout(G.nodes(), shells)
    nx.draw_networkx_nodes(G, pos, node_color=node_colours, node_size=4000, alpha=0.8)
    nx.draw_networkx_edges(G, pos, edge_color=pcolors, edgelist= positive, edge_cmap=plt.cm.Reds, edge_vmin=0.0, width=2.0, alpha=0.4, min_source_margin=25, min_target_margin=25)
    nx.draw_networkx_edges(G, pos, edge_color=ncolors, edgelist= negative, edge_cmap=plt.cm.Blues, edge_vmin=0.0, width=2.0, alpha=0.4, min_source_margin=25, min_target_margin=25)
    nx.draw_networkx_labels(G, pos, font_size=16)
    if knockout == "CTRL":
        plt.title("Spearman correlation network of the control sample", fontsize=16)
    else:
        plt.title(f"Spearman correlation network of the {knockout} knockout experiment", fontsize=16)

    plt.tight_layout()
    fig.savefig(f"../../Images/correlation_network/correlation_network_spearman_{knockout}.pdf")


if __name__ == "__main__":
    systems = pkl.load(open("../../TransformedData_DerivedData/correlation_network/corr_coeff_20.pkl", "rb")).keys()
    with Pool(20) as p:
        p.starmap(corr_to_g, [(s, 0.7, True) for s in systems])