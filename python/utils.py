import numpy as np
import pandas as pd

from os import popen
import ot

Dmax = 200 # constant across comparisons

def sort_dict(d: dict):
    return sorted(d.items(), key=operator.itemgetter(1), reverse=True)

def jaccard_similarity(list_a, list_b):
    set_a = set(list_a)
    set_b = set(list_b)
    numerator = len(set_a.intersection(set_b))
    denominator = len(list_a) + len(list_b) - numerator
    return numerator/denominator

def get_raw_distance_matrix( 
    f1,
    f2,
    as_pandas_dataframe=False,
    index_column=None,
    verbose=True,
    db='/fh/fast/matsen_e/bolson2/transport/iel_data/fake_pubtcrs_db_mouse',
    exe='bin/tcrdists',
    dedup=False
):
    if dedup:
        df_1 = get_df_from_file(f1)
        f1 = "dedup_1.csv"
        df_1.iloc[:, [0, 1]].drop_duplicates().to_csv(f1, header=False, index=False)
        df_2 = get_df_from_file(f2)
        f2 = "dedup_2.csv"
        df_2.iloc[:, [0, 1]].drop_duplicates().to_csv(f2, header=False, index=False)
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

    if 'tcr' not in df.columns:
        df = append_id_column(df)
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
        df = append_id_column(df)
        counter = Counter(df['TCR'])
        unique_tcrs = counter.keys()
        mass = [count/N for count in counter.values()]
        return (mass, unique_tcrs)
    else:
        raise Exception("Unsupported distribution_type")
    return (mass, gene_mass_dict)

def get_transport_objects(
    df_1,
    df_2,
    distribution_type="uniform",
    DMAX=200,
    db='/fh/fast/matsen_e/bolson2/transport/iel_data/fake_pubtcrs_db_mouse',
    exe='bin/tcrdists',
):
    mass_1 = get_mass_objects(df_1, distribution_type=distribution_type)
    mass_2 = get_mass_objects(df_2, distribution_type=distribution_type)

    deduplicated_file_1 = "deduplicated_file_1.csv"
    deduplicated_file_2 = "deduplicated_file_2.csv"

    df_1.iloc[:, [0, 1]].to_csv(deduplicated_file_1, header=False, index=False)
    df_2.iloc[:, [0, 1]].to_csv(deduplicated_file_2, header=False, index=False)

    dist_mat = get_raw_distance_matrix(
        deduplicated_file_1,
        deduplicated_file_2,
        db='/fh/fast/matsen_e/bolson2/transport/iel_data/fake_pubtcrs_db_mouse',
        exe='bin/tcrdists',
        verbose=False)/DMAX

    return mass_1, mass_2, dist_mat

def append_id_column(df):
    if 'TCR' not in df.columns:
        df['TCR'] = [','.join([gene, cdr3]) for gene, cdr3 in zip(df['v_gene'], df['cdr3'])]
    return df

def get_effort_scores(df_1, df_2, LAMBDA=0.1, DMAX=200):
    mass_1, mass_2, dist_mat = get_transport_objects(df_1, df_2)
    ot_mat = ot.sinkhorn(mass_1[0], mass_2[0], dist_mat, LAMBDA)
    effort_mat = np.multiply(dist_mat, ot_mat)

    N2 = df_2.shape[0]
    efforts = DMAX*N2*effort_mat.sum(axis=0)

    assert len(efforts) == N2

    return efforts

