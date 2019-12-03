import glob
import json
import seaborn as sns
import ot

from utils import *

def run_within_gene_analysis(file1, file2, lambd, results_dir, do_clustermap=False, verbose=False):
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
        dist_mat_gene = get_raw_distance_matrix(file1_gene, file2_gene, as_pandas_dataframe=True, index_column=1, verbose=verbose)/Dmax
        if all(dim > 1 for dim in dist_mat_gene.shape):
            ot_mat_gene = pd.DataFrame(ot.sinkhorn(mass_1, 
                                              mass_2, 
                                              dist_mat_gene,
                                              lambd), 
                                  index=df_1_gene.iloc[:, 1], 
                                  columns=df_2_gene.iloc[:, 1]
                                 )
            effort_mat_gene = ot_mat_gene.multiply(dist_mat_gene)
            row_effort_scores = effort_mat_gene.sum(axis=1)
            column_effort_scores = effort_mat_gene.sum(axis=0)
            scores[gene] = {}
            scores[gene]['row'] = row_effort_scores
            scores[gene]['column'] = column_effort_scores

            if do_clustermap:
                cplt = sns.clustermap(effort_mat_gene)
                cplt.savefig(results_dir + str(gene) + 'heatmap.png')
    return {"scores": scores}

def get_within_gene_score_results(file1, file_list, lambd, output_filename):
    scores = []
    for filename in file_list:
        print(filename)
        trial_scores = run_within_gene_analysis(file0, filename, lambd=lambd, results_dir=output_filename + str("_dir/"))['scores']
        max_scores = [np.max(trial_scores[gene]['row']) for gene in trial_scores.keys()]
        print("Max score: " + str(np.max(max_scores)))
        scores.extend(max_scores)

    result = {"max_scores": scores}
    with open(output_filename, 'w') as fp:
        json.dump(result, fp)

    return result



if __name__ == "__main__":
    file_dir = "/fh/fast/matsen_e/bolson2/transport/replicates/"
    file0 = file_dir + "Subject1_aliquot01_reparsed.csv"
    file1 = file_dir + "Subject1_aliquot02_reparsed.csv"
    file2 = file_dir + "Subject3_PBMC_aliquot01_reparsed.csv"

    main_dir = "/home/bolson2/sync/within_gene/"
    #analysis_02 = run_within_gene_analysis(file0, file2, lambd=0.01, results_dir=main_dir + "subject_1_3/")
    #analysis_01 = run_within_gene_analysis(file0, file1, lambd=0.01, results_dir=main_dir + "subject_1_1/")

    file_dir = '/fh/fast/matsen_e/bolson2/transport/replicates/'
    files = [ \
        sorted(glob.glob(file_dir + 'Subject1_*.csv')),
        sorted(glob.glob(file_dir + 'Subject2_*.csv')),
        sorted(glob.glob(file_dir + 'Subject3_*.csv')),
        sorted(glob.glob(file_dir + 'Subject4_*.csv'))
    ]
    files[0].remove(file1)

    results = {i: get_within_gene_score_results(file1, files[i], lambd=0.01, output_filename="within_result_" + str(i)) for i in range(len(files))}
    max_score_dict = {result: results[result]['max_scores'] for result in results}
    with open("within_results.json", 'w') as fp:
        json.dump(results, fp)




    import pdb; pdb.set_trace()
