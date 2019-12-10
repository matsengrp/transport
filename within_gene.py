import glob
import json
import ntpath
import ot
import re
import seaborn as sns

from utils import *

def get_differential_loneliness_scores(effort_mat, self_dist_mat, eps):
    loneliness_scores = {}
    tcrs = effort_mat.index
    for t in tcrs:
        neighbors = []
        for t2 in tcrs:
            dist = self_dist_mat.loc[t, t2]
            try:
                if dist < eps:
                    neighbors.append(t2)
            except ValueError:
                try:
                    if any(dist) < eps:
                        neighbors.append(t2)
                except ValueError:
                    pass

        loneliness_scores[t] = effort_mat.loc[neighbors, :].sum()

    return {"scores": loneliness_scores}

def run_within_gene_analysis(file1, file2, lambd, results_dir, method="marginal", do_clustermap=False, verbose=False):
    #dist_mat = get_raw_distance_matrix(file1, file2, as_pandas_dataframe=True, index_column=1)/Dmax
    df_1 = get_df_from_file(file1)
    df_2 = get_df_from_file(file2)
    shared_genes = set(df_1['v_gene']).intersection(set(df_2['v_gene']))
    cdr3_lengths = []
    scores = {}
    for gene in shared_genes:
        if verbose:
            print(gene)
        df_1_gene = df_1[df_1['v_gene'] == gene]
        file1_gene = "df_1_gene.csv"
        df_1_gene.to_csv(file1_gene, header=False, index=False)
        df_2_gene = df_2[df_2['v_gene'] == gene]
        file2_gene = "df_2_gene.csv"
        df_2_gene.to_csv(file2_gene, header=False, index=False)
        mass_1, gene_mass_dict_1 = get_mass_objects(df_1_gene, "inverse_to_v_gene")
        mass_2, gene_mass_dict_2 = get_mass_objects(df_2_gene, "inverse_to_v_gene")
        dist_mat_gene = get_raw_distance_matrix(file1_gene, file2_gene, as_pandas_dataframe=True, index_column=1, verbose=verbose, db='/fh/fast/matsen_e/bolson2/transport/iel_data/fake_pubtcrs_db_mouse')/Dmax
        if all(dim > 1 for dim in dist_mat_gene.shape):
            ot_mat_gene = pd.DataFrame(ot.sinkhorn(mass_1, 
                                              mass_2, 
                                              dist_mat_gene,
                                              lambd), 
                                  index=df_1_gene.iloc[:, 1], 
                                  columns=df_2_gene.iloc[:, 1]
                                 )
            effort_mat_gene = ot_mat_gene.multiply(dist_mat_gene)
            if method == "marginal":
                row_scores = get_loneliness_scores(effort_mat_gene, margin_index="row", method="ratio") 
                column_scores = get_loneliness_scores(effort_mat_gene, margin_index="column", method="ratio")
            elif method == "differential":
                self_dist_mat_gene = get_raw_distance_matrix(file1_gene, file1_gene, as_pandas_dataframe=True, index_column=1, verbose=verbose)/Dmax
                row_scores = get_differential_loneliness_scores(effort_mat_gene, self_dist_mat_gene, eps=0.25)
            scores[gene] = {}
            scores[gene]['row'] = row_scores
            scores[gene]['column'] = column_scores

            if do_clustermap:
                cplt = sns.clustermap(effort_mat_gene)
                cplt.savefig(results_dir + str(gene) + 'heatmap.png')
    return scores

def get_loneliness_scores(effort_matrix, margin_index, method):
    margins = {"row": 1, "column": 0}
    N1 = effort_matrix.shape[0]
    N2 = effort_matrix.shape[1]
    if method == "ratio":
        scale_ratio = N1/N2 if margin_index == "row" else N2/N1
    elif method == "normalized":
        scale_ratio = 1/(N1*N2)
    else:
        scale_ratio = 1
    scores = [x*scale_ratio for x in list(effort_matrix.sum(axis=margins[margin_index]))]
    return scores

def get_within_gene_score_results(file1, file_list, lambd, output_filename):
    scores = {}
    max_scores = {}
    for filename in file_list:
        file_subject = re.sub('\.tcrs', '', ntpath.basename(filename))
        print(file_subject)
        scores[file_subject] = run_within_gene_analysis(file1, filename, lambd=lambd, results_dir=output_filename + str("_dir/"))

    result = {"scores": scores}
    with open(output_filename, 'w') as fp:
        json.dump(result, fp)

    return result



if __name__ == "__main__":

    data_to_use = "iel"
    if data_to_use == "iel":
        file_dir = '/fh/fast/matsen_e/bolson2/transport/iel_data/iels_tcrs_by_mouse/'
        file1 = file_dir + 'CD4_17_B.tcrs' # Seems to have a moderate number of sequences (455)
        file0 = file_dir + 'CD4_9_B.tcrs'
        file2 = file_dir + 'CD8_10_B.tcrs'
        main_dir = "/home/bolson2/sync/within_gene/"
        #analysis_02 = run_within_gene_analysis(file0, file2, lambd=0.01, results_dir=main_dir + "CD4-CD8/", do_clustermap=True)
        #analysis_01 = run_within_gene_analysis(file0, file1, lambd=0.01, results_dir=main_dir + "CD4-only/", do_clustermap=True)
        files = [ \
            sorted(glob.glob(file_dir + 'CD4*B.tcrs')),
            sorted(glob.glob(file_dir + 'CD8*B.tcrs')),
            sorted(glob.glob(file_dir + 'DN*B.tcrs'))
        ]
        files[2].remove(file1)
    elif data_to_use == "adaptive":
        file_dir = "/fh/fast/matsen_e/bolson2/transport/replicates/"
        file0 = file_dir + "Subject1_aliquot01_reparsed.csv"
        file1 = file_dir + "Subject1_aliquot02_reparsed.csv"
        file2 = file_dir + "Subject3_PBMC_aliquot01_reparsed.csv"

        main_dir = "/home/bolson2/sync/within_gene/"
        file_dir = '/fh/fast/matsen_e/bolson2/transport/replicates/'
        files = [ \
            sorted(glob.glob(file_dir + 'Subject1_*.csv')),
            sorted(glob.glob(file_dir + 'Subject2_*.csv')),
            sorted(glob.glob(file_dir + 'Subject3_*.csv')),
            sorted(glob.glob(file_dir + 'Subject4_*.csv'))
        ]
        files[0].remove(file1)

    results = {i: get_within_gene_score_results(file1, files[i], lambd=0.01, output_filename="within_result_" + str(i)) for i in range(len(files))}
    with open("within_results.json", 'w') as fp:
        json.dump(results, fp)

    import pdb; pdb.set_trace()
