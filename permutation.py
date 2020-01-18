from within_gene import *

def get_trial_scores(file_1, file_2):
    scores = run_within_gene_analysis(file_1, file_2, results_dir=None, method=sys.argv[1])
    return scores

def get_permutation_score_distribution(file_1, file_2, trial_count=50):
    df_1 = get_df_from_file(file_1)
    df_2 = get_df_from_file(file_2)
    N1 = df_1.shape[0]
    N2 = df_2.shape[0]
    full_df = pd.concat([df_1, df_2], axis=0).reset_index(drop=True)
    all_scores = []
    for trial in range(trial_count):
        print(trial)
        df_1_trial_indices = sample(full_df.index.tolist(), N1)
        df_1_trial = full_df.iloc[df_1_trial_indices]
        df_2_trial = full_df.iloc[~full_df.index.isin(df_1_trial_indices)]

        file_1_trial = "df_1_trial.csv"
        file_2_trial = "df_2_trial.csv"
        df_1_trial.to_csv(file_1_trial, header=False, index=False)
        df_2_trial.to_csv(file_2_trial, header=False, index=False)

        trial_scores = get_trial_scores(file_1_trial, file_2_trial)
        all_scores.append(trial_scores)

    return all_scores


if __name__ == "__main__":
    file_dir = '/fh/fast/matsen_e/bolson2/transport/iel_data/iels_tcrs_by_mouse/'
    cd4_df = file_dir + "CD4_17_B.tcrs"
    cd8_df = file_dir + 'CD8_10_B.tcrs'
    dn_df = file_dir + 'DN_15_B.tcrs' 
    
    cd4_cd8_scores = get_permutation_score_distribution(cd4_df, cd8_df)
    cd4_dn_scores = get_permutation_score_distribution(cd4_df, dn_df)
    cd8_dn_scores = get_permutation_score_distribution(cd8_df, dn_df)
    results = [
        {"CD4_CD8": cd4_cd8_scores},
        {"CD4_DN": cd4_dn_scores},
        {"CD8_DN": cd8_dn_scores}
    ]
    
    method = sys.argv[1]
    out_dir = "/home/bolson2/sync/" + method + "/"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    with open(out_dir + "permutation_results.json", 'w') as fp:
        json.dump(results, fp)
