import numpy as np
import ot
import re

from phb_analyses.find_lonely_iels import (
    db,
    Dmax, 
    exe,
    get_raw_distance_matrix,
    lambd,
)

def collapse_allele(gene: str):
    return re.sub("\*[0-9]+", "", gene)

def collapse_gene_subfamily(gene: str):
    return re.sub("\-[0-9]+", "", collapse_allele(gene))

def get_df_from_file(filename, collapse_by_allele=True, collapse_by_subfamily=True):
    df = pd.read_csv(filename, header=None)
    df.columns = ('v_gene', 'cdr3')
    if collapse_by_subfamily:
        df['v_gene'] = [collapse_gene_subfamily(gene) for gene in df['v_gene']]
    elif collapse_by_allele:
        df['v_gene'] = [collapse_allele(gene) for gene in df['v_gene']]
    return df

import pandas as pd
if __name__ == "__main__":
    file1 = "tmptcrs_dir/tmptcrs.0.009219549909406655_f1_tcrs_subset.txt"
    file2 = "tmptcrs_dir/tmptcrs.0.009219549909406655_f2_tcrs_subset.txt"
    df_1 = get_df_from_file(file1)
    df_2 = get_df_from_file(file2)

    N1 = df_1.shape[0]
    N2 = df_2.shape[0]

    mass_1 = np.ones((N1, ))/N1
    mass_2 = np.ones((N2, ))/N2
    dist_mat = get_raw_distance_matrix(file1, file2)/Dmax

    ot_mat = ot.sinkhorn(mass_1, mass_2, dist_mat, lambd)

    gene_transfer_map = {}
    for i in range(0, N1):
        for j in range(0, N2):
            gene_i = df_1.iloc[i, 0]
            gene_j = df_2.iloc[j, 0]
            if gene_i not in gene_transfer_map:
                gene_transfer_map[gene_i] = {}
            if gene_j not in gene_transfer_map[gene_i]:
                gene_transfer_map[gene_i][gene_j] = 0
            # Add the mass transfered from tcr_i to tcr_j to the transfer map for (gene_i, gene_j)
            gene_transfer_map[gene_i][gene_j] += ot_mat[i, j]
