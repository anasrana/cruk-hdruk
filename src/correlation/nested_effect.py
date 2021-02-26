import pandas as pd
import numpy as np
import os
import networkx as nx
import matplotlib.pyplot as plt
import sys

data_dir = '/shared/editter/'

# Get dataframe of knockout genes, ie 20x20
def get_ko_df(input_df):
    return input_df[input_df.index]

# Create binary matrix from 20x20 knockout gene matrix. 
# zero if iloc[i,j] is less than cut-off, else 1
# ie if knockout gene A and gene B counts decrease, assume gene B is also knocked out
def create_gamma_df(input_df, cut_off):
    ko_df = get_ko_df(input_df)
    gamma_df = ko_df.copy(deep=True)
    # Set diagonal elements to 1
    for i in range(0, len(ko_df.index)):
        gamma_df.iloc[i, i] = 1.0
        # 
        for j in range(0, len(ko_df.index)):
            if i != j:
                gamma_df.iloc[i, j] = np.heaviside(-ko_df.iloc[i, j] + cut_off, 1.0)
    return gamma_df

def num_non_zero_entries(df):
    return np.count_nonzero(df)

# Let only the largest DE gene in each column be 1, all others zero
# limitations of this approach
def create_theta_df(input_df):
    theta_df = input_df.copy(deep=True)
    theta_df.iloc[:, :] = 0.0
    for col_idx, col in enumerate(input_df.columns.values):
        idx = np.argmax(np.abs(input_df[col]))
        theta_df.iloc[idx, col_idx] = 1.0
    return theta_df

# 
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

#
def column_filter(input_df, threshold):
    drop_columns = []
    for gene in input_df.columns:
        changed = np.heaviside(np.abs(input_df[gene].values) - threshold, 0.0)
        if not np.any(changed):
            drop_columns.append(gene)
    return drop_columns

# For range of thresholds, give length of array of most expressed genes (for dimensionalty reduction in analysis)
def select_thresh_forGeneNum(num_genes, df):
    # Make dataframe into matrix, then array, sort by descending order. 
    # Values in array will be df values ie here DE 
    df2 = np.abs(df)
    BC = df2.to_numpy()
    BC = np.sort(BC, axis=None)
    BC = BC[::-1]
    thresh_needed = BC[num_genes]
    keep_genes = column_filter_keep(df, thresh_needed)
    return keep_genes

def column_filter_keep(input_df, threshold):
    keep_columns = []
    for gene in input_df.columns:
        changed = np.heaviside(np.abs(input_df[gene].values) - threshold, 0.0)
        if np.any(changed):
            keep_columns.append(gene)
    return keep_columns


def mostDE_genes_for_all_thresh(df):
    threshold_arr = np.geomspace(0.001, 15.0, 300)
    lengths = [len(column_filter(df, thresh)) for thresh in threshold_arr]
    return lengths

def gamma_graph(gamma):
    G = nx.DiGraph()
    G.add_nodes_from(np.arange(0, len(gamma.index), 1))
    labels = {}
    node_colors = []
    colors = []
    sizes = []
    for idx, node in enumerate(G.nodes()):
        labels[node] = gamma.index.values[idx]
        sizes.append(2500)
        node_colors.append('deepskyblue')
    
    for i in range(0, len(gamma.index)):
        for j in range(0, len(gamma.index)):
            if gamma.iloc[i, j] != 0:
                G.add_edge(i, j)
                colors.append(gamma.iloc[i, j])
    return G, labels, colors, sizes, node_colors

def full_graph(F, genes):
    F_df = F[genes]
    G = nx.DiGraph()
    # add nodes from knockout genes
    G.add_nodes_from(np.arange(0, len(F_df.index), 1))
    # add these as other nodes
    G.add_nodes_from(np.arange(len(F_df.index), len(F_df.index) + len(F_df.columns), 1))
    
    labels = {}
    node_colors = []
    colors = []
    sizes = []
    for idx, node in enumerate(G.nodes()):
        # for knockout genes
        if idx < len(F_df.index):
            labels[node] = F_df.index.values[idx]
        # for added genes
        else:
            labels[node] = F_df.columns.values[idx - len(F_df.columns)]
    # For knockout genes        
    for i in range(0, len(F_df.index)):
        sizes.append(2500)
        # node_colors.append(np.sum(knock_df.iloc[:, i]))
        node_colors.append('deepskyblue')
        # add edges between knockout genes (currently removed for clarity)
        
    # i knockout genes, j is columns            
    for i in range(0, len(F_df.index)):
        for j in range(len(F_df.index), len(F_df.index) + len(F_df.columns)):
            G.add_edge(i, j)
            colors.append(F_df.iloc[i, j - len(F_df.columns)])
    for i in range(0, len(F_df.columns)):
        node_colors.append('crimson')
        sizes.append(2500)
    return G, labels, colors, sizes, node_colors

def get_shells(F, genes):
    shells = [list(range(len(F.index), len(F.index) + len(genes), 1)), list(range(0, len(F.index), 1))]
    return shells

def plot_gamma_graph(gamma):
    G, labels, colors, sizes, node_colors = gamma_graph(gamma)
    plt.figure(figsize=(12, 12))
    # have to give for input later, position of nodes
    pos = nx.circular_layout(G)

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, cmap=plt.cm.binary, node_size=sizes, alpha=1.0)
    nx.draw_networkx_edges(G, pos, edge_color=colors, edge_cmap=plt.cm.binary, edge_vmin=0.0, width=2.0, arrowstyle='-|>', alpha=0.4, min_source_margin=25, min_target_margin=25)
    nx.draw_networkx_labels(G, pos, labels, font_size=10)
    plt.savefig('gamma_graph_v2.pdf')
    
def plot_shell_graph(F, genes):
    G, labels, colors, sizes, node_colors = full_graph(F, genes)
    plt.figure(figsize=(20, 20))
    shells = get_shells(F, genes)
    pos = nx.shell_layout(G, shells)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=sizes, alpha=0.8)
    nx.draw_networkx_edges(G, pos, edge_color=colors, edge_cmap=plt.cm.binary, edge_vmin=0.0, width=2.0, arrowstyle='-|>', alpha=0.4, min_source_margin=25, min_target_margin=25)
    nx.draw_networkx_labels(G, pos, labels, font_size=10)
    plt.savefig('shell_graph_v2_{}.pdf'.format(len(genes)))


if __name__ == '__main__':

    start_df = pd.read_csv(data_dir + 'log2FoldChangeData.csv', index_col=0)
    print(start_df.head(19))

    # Selected cut-off to be zero
    gamma = create_gamma_df(start_df, 0.0)
    theta = create_theta_df(start_df)
    F = gamma.dot(theta)
    genes = select_thresh_forGeneNum(20, start_df)
    plot_gamma_graph(gamma)
    plot_shell_graph(F, genes)
    
#     plt.imshow(gamma, cmap='Blues')
#     ax = plt.gca()
#     ax.set_xticks([i for i in range(0, len(gamma.index))])
#     ax.set_yticks([i for i in range(0, len(gamma.index))])
#     ax.set_xticklabels(gamma.index.values)
#     ax.set_yticklabels(gamma.index.values)
#     plt.xticks(rotation=90)
#     plt.savefig('gamma.pdf')
    
#     non_zero_entries = []
#     cut_offs = np.linspace(np.min(np.min(get_ko_df(start_df))), 0.0, 100)
#     for cut_off in cut_offs:
#         non_zero_entries.append(num_non_zero_entries(create_gamma_df(start_df, cut_off)))
#     plt.figure(figsize=(6, 6))
#     plt.plot(cut_offs, non_zero_entries, c='navy', lw=2.0)
#     plt.xlabel('Cut-off')
#     plt.ylabel('# of non-zero entries')
#     plt.savefig('cutoffs.pdf')
