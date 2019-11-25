import glob
import json
import matplotlib.pyplot as plt
import numpy as np
import operator
import ot
import pandas as pd
import random
import re
import seaborn as sns; sns.set()
from os import popen

exe = 'bin/tcrdists'
db = 'data/db'
Dmax = 200 # constant across comparisons
lambd = .01

def sort_dict(d: dict):
    return sorted(d.items(), key=operator.itemgetter(1), reverse=True)

def jaccard_similarity(list_a, list_b):
    set_a = set(list_a)
    set_b = set(list_b)
    numerator = len(set_a.intersection(set_b))
    denominator = len(list_a) + len(list_b) - numerator
    return numerator/denominator

def get_raw_distance_matrix( f1, f2 ):
    cmd = '{} -i {} -j {} -d {} --terse'.format( exe, f1, f2, db )
    print(cmd)
    all_dists = []
    for line in popen(cmd):
        try:
            all_dists.append( [float(x) for x in line.split() ] )
        except ValueError:
            import pdb; pdb.set_trace()
            print(line)
    N1 = len(all_dists)
    N2 = len(all_dists[0])
    for dists in all_dists:
        assert len(dists) == N2

    D = np.array(all_dists)
    print('loaded dists',D.shape)
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

def get_ordered_clustermap(mat, order_by="row"):
    clustermap = sns.clustermap(mat)
    indices = clustermap.dendrogram_row.reordered_ind if order_by == "row" else clustermap.dendrogram_col.reordered_ind
    plt.close('all')
    new_clustermap = sns.heatmap(mat.loc[mat.index[indices], mat.index[indices]], xticklabels=True, yticklabels=True)
    return new_clustermap

def plot_distance_versus_transport(vb_distances, transports, colors, filename, gene_transfer_matrix):
    row_gene_factors = pd.factorize(gene_transfer_matrix.index)
    column_gene_factors = pd.factorize(gene_transfer_matrix.columns)
    cmap = plt.cm.Spectral
    norm = plt.Normalize(
            vmin=np.min(row_gene_factors[0]), 
            vmax=max(row_gene_factors[0])
    )

    plt.subplot(2, 1, 1)
    plt.scatter(vb_distances, transports, c=cmap(norm(pd.factorize(colors)[0])))
    axes = plt.gca()
    axes.set_xlim(-3, 73)
    axes.set_ylim([np.min(transports), np.max(transports)])
    plt.xlabel("Distance between VB genes")
    plt.ylabel("Cumulative transport")
    plt.title("Cumulative transport versus gene distance for lambda = " + str(lambd))

    plt.subplot(2, 1, 2)
    plt.scatter(vb_distances, np.log(transports), c=cmap(norm(pd.factorize(colors)[0])))
    axes = plt.gca()
    axes.set_xlim(-3, 73)
    plt.xlabel("Distance between VB genes")
    plt.ylabel("log(Cumulative transport)")
    plt.title("log(Cumulative transport) versus gene distance for lambda = " + str(lambd))
    
    plt.savefig(filename)

def get_mass_objects(df, distribution_type):
    if distribution_type == "inverse_to_v_gene":
        def get_gene_masses(gene_list):
            unique_genes = list(set(gene_list))
            gene_mass_dict = {gene: 1/len(unique_genes) for gene in unique_genes}
            return gene_mass_dict

        mass = get_gene_weighted_mass_distribution(df)
        gene_mass_dict = get_gene_masses(df['v_gene'])
    elif distribution_type == "uniform":
        N = df.shape[0]
        mass = np.ones((N, ))/N
        gene_mass_dict = tabulate_gene_frequencies(list(df['v_gene']))
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

def run_gene_score_analysis(file1, file2, do_plot=False):
    df_1 = get_df_from_file(file1)
    df_2 = get_df_from_file(file2)

    N1 = df_1.shape[0]
    N2 = df_2.shape[0]

    results = {}
    score_matrices = {}
    gene_transfer_matrices = {}
    distribution_types = ["inverse_to_v_gene"]
    lambdas = [0.01]
    for distribution_type in distribution_types:
        mass_1, gene_mass_dict_1 = get_mass_objects(df_1, distribution_type)
        mass_2, gene_mass_dict_2 = get_mass_objects(df_2, distribution_type)

        dist_mat = get_raw_distance_matrix(file1, file2)/Dmax

        for lambd in lambdas:
            results[lambd] = {}

            ot_mat = ot.sinkhorn(mass_1, mass_2, dist_mat, lambd)

            gene_transfer_matrix = get_gene_transfer_matrix(df_1, df_2, ot_mat)
            gene_transfer_matrices[lambd] = gene_transfer_matrix


            # Get distance matrix corresponding entrywise to gene_transfer_matrix 
            vb_mat_sub = vb_mat.loc[gene_transfer_matrix.index, gene_transfer_matrix.columns]

            (vb_distances, transports, column_colors, score_matrix) = get_scoring_quantities(gene_transfer_matrix, vb_mat_sub)

            results[lambd]['score_matrix'] = score_matrix
            results[lambd]['scores'] = dict(score_matrix.sum(axis=0))
            results[lambd]['jaccard_index'] = jaccard_similarity(gene_transfer_matrix.index, gene_transfer_matrix.columns)

            score_matrices[lambd] = score_matrix

        if do_plot:
            plot_distance_versus_transport(vb_distances, transports, column_colors, results_dir + "scatterplot_grid.png", gene_transfer_matrix)

    for lambd in lambdas:
        score_plt = sns.clustermap(score_matrices[lambd], xticklabels=True, yticklabels=True)
        score_plt.savefig(results_dir + "score_" + str(lambd) + ".png")

    for lambd in lambdas:
        gene_transfer_matrix = gene_transfer_matrices[lambd]
        transfer_plt = get_ordered_clustermap(gene_transfer_matrix)
        fig = transfer_plt.get_figure()
        fig.savefig(results_dir + "gene_transfer_" + str(lambd) + ".png")
        plt.close()
        #gene_transfer_matrix = gene_transfer_matrix[gene_transfer_matrix.index]
        transfer_plt = sns.clustermap(gene_transfer_matrix, xticklabels=True, yticklabels=True)
        transfer_plt.savefig(results_dir + distribution_type + "_transfer_heatmap_" + "lambda_" + str(lambd) + ".png")


    return results

def get_score_results(file1, file_list, output_filename):
    max_scores = []
    jaccard_indices = []
    for filename in file_list:
        print(filename)
        trial_scores = run_gene_score_analysis(file1, filename, do_plot=False)[0.01]
        max_score = np.max(list(trial_scores["scores"].values()))
        print("Max score: " + str(max_score) + ", Jaccard index: " + str(trial_scores['jaccard_index']))
        max_scores.append(max_score)
        jaccard_indices.append(trial_scores['jaccard_index'])
    result = {"max_scores": max_scores, "jaccard_indices": jaccard_indices}
    with open(output_filename, 'w') as fp:
        json.dump(result, fp)
    return result

if __name__ == "__main__":
    va_mat = get_gene_distance_matrix("data/gene_dist_matrices/va_dist.txt")
    vb_mat = get_gene_distance_matrix("data/gene_dist_matrices/vb_dist.txt")

    results_dir = "results/gene_transfer/"
    do_full = False
    use_replicate_data = True
    if do_full:
        file1 = "data/iel_data/ielrep_beta_CD4_tcrs.txt" 
        file2 ="data/iel_data/ielrep_beta_DN_tcrs.txt"  
    elif use_replicate_data:
        file_dir = "/fh/fast/matsen_e/bolson2/transport/replicates/"
        file0 = file_dir + "Subject1_aliquot02_reparsed.csv" 
        file1 = file_dir + "Subject1_aliquot01_reparsed.csv" 
        file2 = file_dir + "Subject2_aliquot01_reparsed.csv" 
        file3 = file_dir + "Subject3_PBMC_aliquot01_reparsed.csv" 
        file4 = file_dir + "Subject4_PBMC_aliquot01_reparsed.csv" 
    else:
        file1 = "data/yfv/P1_0_F1_.txt.top1000.tcrs"
        file2 = "data/yfv/P1_0_F2_.txt.top1000.tcrs"


    if True:
        obs_scores_01 = run_gene_score_analysis(file1, file0, do_plot=True)[0.01]['scores']
        obs_scores_12 = run_gene_score_analysis(file1, file2, do_plot=False)[0.01]['scores']
        obs_scores_13 = run_gene_score_analysis(file1, file3, do_plot=False)[0.01]['scores']
        obs_scores_14 = run_gene_score_analysis(file1, file4, do_plot=False)[0.01]['scores']
        obs_scores_12_sorted = sort_dict(obs_scores_12)
        obs_scores_13_sorted = sort_dict(obs_scores_13)
        obs_scores_14_sorted = sort_dict(obs_scores_14)
        print(obs_scores_12_sorted[0]) 
        print(obs_scores_13_sorted[0])
        print(obs_scores_14_sorted[0])

    file_dir = '/fh/fast/matsen_e/bolson2/transport/replicates/'
    files = [ \
        sorted(glob.glob(file_dir + 'Subject1_*.csv')),
        sorted(glob.glob(file_dir + 'Subject2_*.csv')),
        sorted(glob.glob(file_dir + 'Subject3_*.csv')),
        sorted(glob.glob(file_dir + 'Subject4_*.csv'))
    ]
    files[0].remove(file1)

    results = {i: get_score_results(file1, files[i], "result_boot_" + str(i)) for i in range(len(files))}
    max_score_dict = {result: results[result]['max_scores'] for result in results}
    jaccard_index_dict = {result: results[result]['jaccard_indices'] for result in results}
    with open("results_data_boot.json", 'w') as fp:
        json.dump(results, fp)
    import pdb; pdb.set_trace()

