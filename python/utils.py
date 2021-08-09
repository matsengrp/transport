import numpy as np
import pandas as pd

import os
import ot

from python.tcr_dist import TCRDist
from config import CONFIG

def sort_dict(d: dict):
    return sorted(d.items(), key=operator.itemgetter(1), reverse=True)

def jaccard_similarity(list_a, list_b):
    set_a = set(list_a)
    set_b = set(list_b)
    numerator = len(set_a.intersection(set_b))
    denominator = len(list_a) + len(list_b) - numerator
    return numerator/denominator

def collapse_allele(gene: str):
    return re.sub("\*[0-9]+", "", gene)

def collapse_gene_subfamily(gene: str):
    return re.sub("\-[0-9]+", "", collapse_allele(gene))

def get_df_from_file(filename, collapse_by_allele=False, collapse_by_subfamily=False):
    df = pd.read_csv(filename, header=None)
    df.columns = ['v_gene', 'cdr3']
    if collapse_by_subfamily:
        df['v_gene'] = [collapse_gene_subfamily(gene) for gene in df['v_gene']]
    elif collapse_by_allele:
        df['v_gene'] = [collapse_allele(gene) for gene in df['v_gene']]

    if 'tcr' not in df.columns:
        df = append_id_column(df)
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

def collapse_duplicates(v):
    starting_indices = [0]
    for i in range(1, len(v)):
        if v[i] != v[i - 1]:
            starting_indices.append(i)
    return starting_indices

def write_deduplicated_file(df, filename, output_dir=CONFIG["TMP_OUTPUT"]):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df.iloc[:, [0, 1]].drop_duplicates().to_csv(os.path.join(output_dir, filename), header=False, index=False)


def append_id_column(df):
    if 'tcr' not in df.columns:
        df['tcr'] = [','.join([gene, cdr3]) for gene, cdr3 in zip(df['v_gene'], df['cdr3'])]
    return df

def extract_cdr3s(sequences):
    return [s.split(',')[1] for s in sequences]
