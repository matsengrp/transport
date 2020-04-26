import glob
import json
import ntpath
import os
import ot
import re
import seaborn as sns
import sys

from random import sample

from utils import *

LAMBDA = 0.1

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

def downsample_dfs(df_1, df_2):
    N1 = df_1.shape[0]
    N2 = df_2.shape[0]
    if N1 > N2:
        df_1 = df_1.loc[sample(df_1.index.tolist(), N2)]
    else:
        df_2 = df_2.loc[sample(df_2.index.tolist(), N1)]
    assert df_1.shape[0] == df_2.shape[0]
    
    return df_1, df_2

def run_within_gene_analysis(file1, file2, results_dir, method, lambd=LAMBDA, do_clustermap=False, verbose=False):
    df_1 = get_df_from_file(file1)
    df_2 = get_df_from_file(file2)

    shared_genes = set(df_1['v_gene']).intersection(set(df_2['v_gene']))
    cdr3_lengths = []
    scores = {}
    for gene in shared_genes:
        if verbose:
            print(gene)
        df_1_gene = df_1[df_1['v_gene'] == gene]
        df_2_gene = df_2[df_2['v_gene'] == gene]
        if method in ("subsample", "subsample_avg"):
            df_1_gene, df_2_gene = downsample_dfs(df_1_gene, df_2_gene)

        file1_gene = "df_1_gene.csv"
        file2_gene = "df_2_gene.csv"
        df_1_gene.to_csv(file1_gene, header=False, index=False)
        df_2_gene.to_csv(file2_gene, header=False, index=False)
        mass_1, gene_mass_dict_1 = get_mass_objects(df_1_gene, "inverse_to_v_gene")
        mass_2, gene_mass_dict_2 = get_mass_objects(df_2_gene, "inverse_to_v_gene")
        dist_mat_gene = get_raw_distance_matrix(file1_gene, file2_gene, as_pandas_dataframe=True, index_column=1, verbose=verbose, db='/fh/fast/matsen_e/bolson2/transport/iel_data/fake_pubtcrs_db_mouse')/Dmax
        if all(dim > 5 for dim in dist_mat_gene.shape):
            ot_mat_gene = pd.DataFrame(ot.sinkhorn(mass_1, 
                                              mass_2, 
                                              dist_mat_gene,
                                              lambd), 
                                  index=df_1_gene.iloc[:, 1], 
                                  columns=df_2_gene.iloc[:, 1]
                                 )
            effort_mat_gene = ot_mat_gene.multiply(dist_mat_gene)
            if method == "differential":
                self_dist_mat_gene = get_raw_distance_matrix(file1_gene, file1_gene, as_pandas_dataframe=True, index_column=1, verbose=verbose)/Dmax
                row_scores = get_differential_loneliness_scores(effort_mat_gene, self_dist_mat_gene, eps=0.25)
            else:
                row_scores = get_loneliness_scores(effort_mat_gene, margin_index="row", method=method) 
                column_scores = get_loneliness_scores(effort_mat_gene, margin_index="column", method=method)
            scores[gene] = {}
            scores[gene]['row'] = row_scores
            scores[gene]['column'] = column_scores

            if do_clustermap:
                cplt = sns.clustermap(effort_mat_gene)
                cplt.savefig(results_dir + str(gene) + 'heatmap.png')
    return scores

def get_loneliness_scores(effort_matrix, margin_index, method):
    margins = {"row": 1, "column": 0}
    cdr3s = {"row": effort_matrix.index, "column": effort_matrix.columns}
    N1 = effort_matrix.shape[0]
    N2 = effort_matrix.shape[1]
    if method == "ratio":
        scale_ratio = N1 if margin_index == "row" else N2
    elif method == "normalized":
        scale_ratio = 1/(N1*N2)
    elif method == "subsample":
        scale_ratio = N1
    else:
        scale_ratio = 1
    if method == "subsample_avg":
        scores = [x*scale_ratio for x in list(effort_matrix.mean(axis=margins[margin_index]))]
    else:
        effort_sums = effort_matrix.sum(axis=margins[margin_index])
        scores = list(effort_sums*scale_ratio) 
        cdr3s = effort_sums.index
        #scores = {cdr3: (effort_sums[cdr3]).item()*scale_ratio for cdr3 in effort_sums}
        try:
            scores = [{cdr3s[i]: scores[i]} for i in range(len(effort_sums))]
        except:
            import pdb; pdb.set_trace()
        #scores = {cdr3: scale_ratio*effort_sum.item() for cdr3, effort_sum in effort_sums.items()}
        #scores = {x*scale_ratio for x in list(effort_matrix.sum(axis=margins[margin_index]))]
    return scores

def get_within_gene_score_results(file1, file_list, method, lambd, output_filename):
    scores = {}
    for file_1 in file1:
        for file_2 in file_list:
            if file_1 is not file_2:
                file_subject = re.sub('\.tcrs', '', ntpath.basename(file_2))
                print(file_subject)
                scores[file_subject] = run_within_gene_analysis(file_1, file_2, lambd=lambd, results_dir=output_filename + str("_dir/"), method=method)
            else:
                print("Skipping " + file_1 + "...")

        result = [{"biological": {0: scores}}]
    with open(output_filename, 'w') as fp:
        json.dump(result, fp)

    return result



if __name__ == "__main__":

    data_to_use = "iel"
    if data_to_use == "iel":
        file_dir = '/fh/fast/matsen_e/bolson2/transport/iel_data/iels_tcrs_by_mouse/'
        reference_file = sys.argv[1]
        file1 = file_dir + reference_file + '_B.tcrs'
        main_dir = "/home/bolson2/sync/within_gene/"
        files = [ \
            sorted(glob.glob(file_dir + 'CD4*B.tcrs')),
            sorted(glob.glob(file_dir + 'CD8*B.tcrs')),
            sorted(glob.glob(file_dir + 'DN*B.tcrs')),
        ]

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

    method = sys.argv[2]
    reference_group = re.sub("_.*", "", reference_file)
    out_dir = method + "/" + reference_group
    if not os.path.exists(method):
        os.mkdir(method) 
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    cell_groups = ("DN")
    results = {cell_groups[i]: get_within_gene_score_results(files[2], files[i], method=method, lambd=LAMBDA, output_filename="within_result_" + str(i)) for i in range(len(files))}
    with open(out_dir + "/within_results.json", 'w') as fp:
        json.dump(results, fp)

