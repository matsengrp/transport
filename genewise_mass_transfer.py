import numpy as np
import ot
import pandas as pd
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

def get_df_from_file(filename, collapse_by_allele=True, collapse_by_subfamily=False):
    df = pd.read_csv(filename, header=None)
    df.columns = ('v_gene', 'cdr3')
    if collapse_by_subfamily:
        df['v_gene'] = [collapse_gene_subfamily(gene) for gene in df['v_gene']]
    elif collapse_by_allele:
        df['v_gene'] = [collapse_allele(gene) for gene in df['v_gene']]
    return df

def tabulate_gene_frequencies(gene_list, as_probability=True):
    gene_list_len = len(gene_list)
    unique_genes = list(set(gene_list))
    if as_probability:
        frequency_dict = {gene: gene_list.count(gene)/gene_list_len for gene in unique_genes}
    else:
        frequency_dict = {gene: gene_list.count(gene) for gene in unique_genes}
    return frequency_dict

def get_gene_weighted_mass_distribution(df):
    gene_freqs =  tabulate_gene_frequencies(list(df['v_gene']), as_probability=False)
    num_genes = len(gene_freqs)
    mass_distribution = [1/(num_genes*gene_freqs[gene]) for gene in df['v_gene']]
    return mass_distribution

if __name__ == "__main__":
    file1 = "data/iel_data/ielrep_beta_CD4_tcrs.txt" #"tmptcrs_dir/tmptcrs.0.009219549909406655_f1_tcrs_subset.txt"
    file2 ="data/iel_data/ielrep_beta_DN_tcrs.txt"  #"tmptcrs_dir/tmptcrs.0.009219549909406655_f2_tcrs_subset.txt"
    df_1 = get_df_from_file(file1)
    df_2 = get_df_from_file(file2)

    N1 = df_1.shape[0]
    N2 = df_2.shape[0]

    gene_mass_dict_1 = tabulate_gene_frequencies(list(df_1['v_gene']))
    gene_mass_dict_2 = tabulate_gene_frequencies(list(df_2['v_gene']))

    
    weight_by_v_genes = True
    if weight_by_v_genes:
        mass_1 = get_gene_weighted_mass_distribution(df_1)
        mass_2 = get_gene_weighted_mass_distribution(df_2)
    else:
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
    scores = {}
    scores_min = {}
    for gene in gene_transfer_map:
        scores[gene] = gene_transfer_map[gene][gene]/(gene_mass_dict_1[gene]*gene_mass_dict_2[gene])
        scores_min[gene] = gene_transfer_map[gene][gene]/min(gene_mass_dict_1[gene], gene_mass_dict_2[gene])
    import pdb; pdb.set_trace()
