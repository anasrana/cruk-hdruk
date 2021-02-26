# construct the networks and calculate betweenness centrality
# based on spearman correlation and a 0.7 cutoff

import pandas as pd
import pickle as pkl
import numpy as np
import networkx as nx
from multiprocessing import Manager, Pool

def construct_networks_centrality(system, cut):
    df = pkl.load(open(f"../../TransformedData_DerivedData/correlation_network/corr_coeff_spearman_{system}_KO.pkl", "rb"))
    genes = list(df.columns)
    G = nx.Graph()
    for x, y in np.argwhere(abs(np.triu(df.values)) > cut):
        if x!=y:
            G.add_edge(genes[x], genes[y], weight=df.loc[genes[x], genes[y]])
    pkl.dump(G, open(f"../../TransformedData_DerivedData/correlation_network/graph_spearman_{system}.pkl", "wb"))
    pkl.dump(nx.betweenness_centrality(G, weight=None), open(f"../../TransformedData_DerivedData/correlation_network/betweenness_centrality_spearman_{system}.pkl", "wb"))


if __name__ == "__main__":
    systems = pkl.load(open("../../TransformedData_DerivedData/correlation_network/corr_coeff_20.pkl", "rb")).keys()
    graphs = {}
    bcs = {}
    with Pool(20) as p:
        p.starmap(construct_networks_centrality, [(s, 0.7) for s in systems])


