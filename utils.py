import numpy as np
import pandas as pd

from os import popen

db = '/fh/fast/matsen_e/bolson2/transport/iel_data/fake_pubtcrs_db_mouse'
exe = 'bin/tcrdists'
Dmax = 200 # constant across comparisons

def sort_dict(d: dict):
    return sorted(d.items(), key=operator.itemgetter(1), reverse=True)

def jaccard_similarity(list_a, list_b):
    set_a = set(list_a)
    set_b = set(list_b)
    numerator = len(set_a.intersection(set_b))
    denominator = len(list_a) + len(list_b) - numerator
    return numerator/denominator

def get_raw_distance_matrix( f1, f2, as_pandas_dataframe=False, index_column=None, verbose=True, db=db, exe=exe):
    cmd = '{} -i {} -j {} -d {} --terse'.format( exe, f1, f2, db )
    if verbose:
        print(cmd)
    all_dists = []
    for line in popen(cmd):
        try:
            all_dists.append( [float(x) for x in line.split() ] )
        except ValueError:
            print(line)
    N1 = len(all_dists)
    N2 = len(all_dists[0])
    for dists in all_dists:
        assert len(dists) == N2

    D = np.array(all_dists)
    if verbose:
        print('loaded dists',D.shape)

    if as_pandas_dataframe:
        if index_column:
            D = pd.DataFrame(D)

            df_1 = get_df_from_file(f1)
            df_2 = get_df_from_file(f2)

            import pdb; pdb.set_trace()
            D.index = df_1.iloc[:, index_column]
            D.columns = df_2.iloc[:, index_column]

    return D

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
    return df.drop_duplicates()

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

def get_gene_distance_matrix(filename, collapse_alleles=False):
    df = pd.read_csv(filename, header=None, sep=" ")
    gene_names = df.iloc[:, 1].values
    df = df.drop(columns=[0, 1])

    if collapse_alleles:
        gene_names = [collapse_allele(gene_name[0]) for gene_name in gene_names]
        indices = collapse_duplicates(gene_names)
        df = df.iloc[indices, indices]
        df.columns = [gene_names[i] for i in indices]
    else:
        df.columns = gene_names

    df.index = df.columns
    return df

def get_mass_objects(df, distribution_type):
    if distribution_type == "inverse_to_v_gene":
        def get_gene_masses(gene_list):
            unique_genes = list(set(gene_list))
            gene_mass_dict = {gene: 1/len(unique_genes) for gene in unique_genes}
            return gene_mass_dict

        mass = get_gene_weighted_mass_distribution(df)
        gene_mass_dict = get_gene_masses(df['v_gene'])
        return (mass, gene_mass_dict)
    elif distribution_type == "uniform":
        from collections import Counter
        N = df.shape[0]
        counter = Counter(df['TCR'])
        unique_tcrs = counter.keys()
        mass = [count/N for count in counter.values()]
        return (mass, unique_tcrs)
    else:
        raise Exception("Unsupported distribution_type")
    return (mass, gene_mass_dict)

def get_gene_transfer_matrix(df_1, df_2, ot_mat):
    N1 = df_1.shape[0]
    N2 = df_2.shape[0]
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

    gene_transfer_matrix = pd.DataFrame(gene_transfer_map)
    return(gene_transfer_matrix)

def get_scoring_quantities(gene_transfer_matrix, vb_matrix):
    vb_distances = []
    transports = []
    row_colors = []
    column_colors = []
    scores = {}
    for i in range(0, vb_matrix.shape[0]):
        row_gene = vb_matrix.index[i]
        scores[row_gene] = {}
        for j in range(0, vb_matrix.shape[1]):
            column_gene = vb_matrix.columns[j]
            vb_dist_ij = vb_matrix.iloc[i, j]
            transport_ij = gene_transfer_matrix.iloc[i, j]
            vb_distances.append(vb_dist_ij)
            transports.append(transport_ij)
            row_colors.append(row_gene)
            column_colors.append(column_gene)
            scores[row_gene][column_gene] = transport_ij*vb_dist_ij

    score_matrix = pd.DataFrame(scores)
    return (vb_distances, transports, column_colors, score_matrix)
