import pandas as pd
import argparse
import numpy as np
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler

parser = argparse.ArgumentParser()

parser.add_argument('-gene_list',
                    dest='gene_list',
                    help='List of the used genes')

parser.add_argument('-connectivity_list',
                    dest='conn_list',
                    help='CSV file that describes how the different genes \
                    interact with each other')

parser.add_argument('-data',
                    dest='data',
                    help='RNA-seq data')
parser.add_argument('-save_dir',
                    dest='save_dir',
                    help='Path to save the histograms')


args = parser.parse_args()

if __name__ == '__main__':
    df_conn = pd.read_csv(args.conn_list, index_col=0)
    df_data = pd.read_csv(args.data, index_col=0)
    df_data.set_index('gene', inplace=True)

    # Get the data from Gene A
    df_gene_a = df_data.loc[df_conn['GeneA']]
    df_gene_b = df_data.loc[df_conn['GeneB']]

    # Because of memory issues I am breaking this loops into two iterations and saving two files for all the images
    iterations = 4
    iter_len = len(df_conn) / iterations
    for iter_x in range(iterations):
        histograms = []
        idx_low = int(iter_x * iter_len)
        idx_high = int((iter_x + 1) * iter_len)
        print(f'Currently running analysis{iter_x}')
        print(f'idx: {idx_low}, {idx_high}')
        for x in tqdm(range(idx_low, idx_high)):
            H, _, _ = np.histogram2d(df_gene_a.iloc[x], df_gene_b.iloc[x], bins=32)
            # normalise histogram
            scaler = StandardScaler()
            hist_scaled = scaler.fit_transform(H)
            histograms.append(hist_scaled)
        # Transfom the data into the format that will be used to train the NN
        xx = np.array(histograms)[:, :, :, np.newaxis]
        np.save(args.save_dir + f'/Nxdata_tf_{iter_x}.npy', xx)
        del xx, histograms
