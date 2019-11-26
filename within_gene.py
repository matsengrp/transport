import seaborn as sns
import ot

from utils import *

def run_within_gene_analysis(file1, file2, lambd, results_dir):
    dist_mat = get_raw_distance_matrix(file1, file2, as_pandas_dataframe=True, index_column=1)/Dmax
    df_1 = get_df_from_file(file1)
    df_2 = get_df_from_file(file2)
    mass_1, gene_mass_dict_1 = get_mass_objects(df_1, "inverse_to_v_gene")
    mass_2, gene_mass_dict_2 = get_mass_objects(df_2, "inverse_to_v_gene")
    ot_mat = pd.DataFrame(ot.sinkhorn(mass_1, mass_2, dist_mat, lambd), 
                                      index=df_1.iloc[:, 1], 
                                      columns=df_2.iloc[:, 1]
                                     )
    effort_mat = ot_mat.multiply(dist_mat)
    shared_genes = set(df_1['v_gene']).intersection(set(df_2['v_gene']))
    for gene in shared_genes:
        effort_mat_gene = effort_mat.iloc[df_1[df_1['v_gene'] == gene].index, 
                              df_2[df_2['v_gene'] == gene].index]
        if all(dim > 1 for dim in effort_mat_gene.shape):
            cplt = sns.clustermap(effort_mat_gene)
            cplt.savefig(results_dir + str(gene) + 'heatmap.png')


if __name__ == "__main__":
    file_dir = "/fh/fast/matsen_e/bolson2/transport/replicates/"
    file0 = file_dir + "Subject1_aliquot01_reparsed.csv"
    file1 = file_dir + "Subject1_aliquot02_reparsed.csv"
    file2 = file_dir + "Subject3_PBMC_aliquot01_reparsed.csv"

    main_dir = "/home/bolson2/sync/within_gene/"
    analysis_02 = run_within_gene_analysis(file0, file2, lambd=0.01, results_dir=main_dir + "subject_1_3/")
    analysis_01 = run_within_gene_analysis(file0, file1, lambd=0.01, results_dir=main_dir + "subject_1_1/")


    import pdb; pdb.set_trace()
