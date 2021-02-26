def create_ko_graph(input_df):
    # selects knockout genes from dataframe, get 20x20df
    knock_df = np.abs(knock_out_df(input_df))
    # reate directional graph (with arrows on edge lines)
    G = nx.DiGraph()
    # nodes labeled by integers. will become object to refer to later
    G.add_nodes_from(np.arange(0, len(knock_df.columns), 1))
    
    #labels=gene name
    labels = {}
    #colour of node
    node_colors = []
    #colour of edge line, high = dark
    colors = []
    sizes = []
    
    # Assign labels (gene names) to node (through dictionary)
    for idx, node in enumerate(G.nodes()):
        labels[node] = knock_df.columns.values[idx]
    # Finding size of Node (bigger if more DE genes, sum across row)
    # square matrix, use i as row and column
    for i in range(0, len(knock_df.columns)):
        sizes.append(1000 * np.sum(knock_df.iloc[i, :]))
        # Node colour is sum of columns
        node_colors.append(np.sum(knock_df.iloc[:, i]))
        #looping over columns
        for j in range(0, len(knock_df.columns)):
            # ifnot diagonal element
            if i != j:
                G.add_edge(i, j)
                # colour of line is abs value of DE
                colors.append(knock_df.iloc[i, j])
    return G, labels, colors, sizes, node_colors

def full_graph(input_df, num):
    knock_df = np.abs(knock_out_df(input_df))
    G = nx.DiGraph()
    # add nodes from knockout genes
    G.add_nodes_from(np.arange(0, len(knock_df.columns), 1))
    # Add extra genes from select_thresh_forGeneNum
    columns = select_thresh_forGeneNum(num, input_df)
    add_df = input_df[columns]
    # add these as other nodes
    G.add_nodes_from(np.arange(len(knock_df.columns), len(add_df.index) + len(add_df.columns), 1))
    
    labels = {}
    node_colors = []
    colors = []
    sizes = []
    for idx, node in enumerate(G.nodes()):
        # for knockout genes
        if idx < len(knock_df.columns):
            labels[node] = knock_df.columns.values[idx]
        # for added genes
        else:
            labels[node] = add_df.columns.values[idx - len(knock_df.columns)]
    # For knockout genes        
    for i in range(0, len(knock_df.columns)):
        sizes.append(2500)
        # node_colors.append(np.sum(knock_df.iloc[:, i]))
        node_colors.append('forestgreen')
        # add edges between knockout genes (currently removed for clarity)
        for j in range(0, len(knock_df.columns)):
            if i != j:
                continue
                #G.add_edge(i, j)
                #colors.append(knock_df.iloc[i, j])
    # i knockout genes, j is columns            
    for i in range(0, len(knock_df.columns)):
        for j in range(len(knock_df.columns), len(knock_df.columns) + len(add_df.columns)):
            G.add_edge(i, j)
            colors.append(add_df.iloc[i, j - len(knock_df.columns)])
    for i in range(0, len(add_df.columns)):
        node_colors.append('indianred')
        sizes.append(2500)
    return G, labels, colors, sizes, node_colors

def get_shells(input_df, num):
    columns = select_thresh_forGeneNum(num, input_df)
    shells = [list(range(len(input_df.index), len(input_df.index) + len(columns), 1)), list(range(0, len(input_df.index), 1))]
    return shells

def plot_ko_graph(input_df):
    G, labels, colors, sizes, node_colors = create_ko_graph(start_df)
    plt.figure(figsize=(12, 12))
    # have to give for input later, position of nodes
    pos = nx.circular_layout(G)

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, cmap=plt.cm.Greens, node_size=sizes, alpha=1.0)
    nx.draw_networkx_edges(G, pos, edge_color=colors, edge_cmap=plt.cm.binary, edge_vmin=0.0, width=2.0, arrowstyle='-|>', alpha=0.4, min_source_margin=25, min_target_margin=25)
    nx.draw_networkx_labels(G, pos, labels, font_size=10)
    plt.savefig('/shared/ko_graph.pdf')
    
def plot_shell_graph(input_df, num):
    G, labels, colors, sizes, node_colors = full_graph(start_df, num)
    plt.figure(figsize=(20, 20))
    shells = get_shells(start_df, num)
    pos = nx.shell_layout(G, shells)
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=sizes, alpha=0.8)
    nx.draw_networkx_edges(G, pos, edge_color=colors, edge_cmap=plt.cm.binary, edge_vmin=0.0, width=2.0, arrowstyle='-|>', alpha=0.4, min_source_margin=25, min_target_margin=25)
    nx.draw_networkx_labels(G, pos, labels, font_size=10)
    plt.savefig('/shared/shell_graph_{}.pdf'.format(num))
    

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

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

# Generate df for knock out genes
def knock_out_df(input_df):
    return input_df[input_df.index]

# Defining function to compare control to mcf7 data, binarise above threshold to 1
def column_filter(input_df, threshold):
    drop_columns = []
    for gene in input_df.columns:
        changed = np.heaviside(np.abs(input_df[gene].values) - threshold, 0.0)
        if not np.any(changed):
            drop_columns.append(gene)
    return drop_columns

def df_to_network_binary(input_df, threshold):
    control = input_df['CTRL'].values
    df = input_df.drop(['CTRL'], axis=1)
    for gene in df.columns:
        df[gene] = np.heaviside(np.abs(df[gene] - control) - threshold, 0.0)
    return df

# Defining function to compare control to mcf7 data, binarise above threshold to 1
def df_to_network_binary2(input_df, threshold):
    for gene in df.columns:
        df[gene] = np.heaviside(np.abs(df[gene]) - threshold, 0.0)
    return df

# Defining function to compare control to mcf7 data, binarise above threshold to 1
def df_to_network_difference(input_df):
    control = input_df['CTRL'].values
    df = input_df.drop(['CTRL'], axis=1)
    for gene in df.columns:
        df[gene] = df[gene] - control
    return df

def plot_heatmap(df, flag_heatmap):
    if flag_heatmap==1:
        fig = plt.figure(figsize=(8, 8))
        sns.heatmap(df)
        plt.show()

def make_Net(df):
    G = nx.Graph()
    G.add_nodes_from(np.arange(0, len(df.columns), 1))
    for i in range(0, len(df.columns)):
        for j in range(0, len(df[df.columns[0]].values)):
            if df.iloc[i, j] == 1.0:
                G.add_edge(i, j)
    return G        
        
def plot_NetGraph(G, df, flag):
    if flag==1:
        fig = plt.figure(figsize=(10, 10))
        labels = {}
        for idx, node in enumerate(G.nodes()):
            labels[node] = df.columns.values[idx]
        pos = nx.circular_layout(G)
        plt.figure(figsize=(10, 10))
        nx.draw_networkx_nodes(G, pos,
                            node_color='#419D78',
                            node_size=2500,
                            alpha=1.0)
        edges = nx.draw_networkx_edges(G, pos,
                            edge_color='navy',
                            width=1.5, 
                            alpha=0.75)
        nx.draw_networkx_labels(G, pos, labels, font_size=12)
        plt.show()

def make_all(input_df, threshold, flag_heatmap, flag_netGraph):
    df = df_to_network_binary(input_df, threshold)
    plot_heatmap(df, flag_heatmap)
    G = make_Net(df)
    plot_NetGraph(G, df, flag_netGraph)
