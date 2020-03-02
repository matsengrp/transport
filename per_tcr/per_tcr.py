from collections import defaultdict
import json
import os
from random import sample
import sys

import ot

sys.path.append('..')
from utils import *

LAMBDA = 0.01
DMAX = 200

def append_id_column(df, prefix):
    df['id'] = ['_'.join([prefix, str(i)]) for i in df.index]
    return df

def get_processed_df(subject, file_dir):
    filename = file_dir + subject + '_B.tcrs'
    df = get_df_from_file(filename)
    df = append_id_column(df, subject)
    return df

def split_datasets(full_df, N1, N2, tcr_index):
    indices_to_split = full_df.index.tolist()
    indices_to_split.remove(tcr_index) # Ensure tcr_index is not in df_1_trial, and necessarily is in df_2_trial
    df_1_trial_indices = sample(indices_to_split, N1)
    df_1_trial = full_df.iloc[df_1_trial_indices]
    df_2_trial = full_df.iloc[~full_df.index.isin(df_1_trial_indices)]
    return df_1_trial, df_2_trial

def do_randomization_test(df_1, df_2, method_type, trial_count=2):
    full_df = pd.concat([df_1, df_2], axis=0).reset_index(drop=True)
    N1 = df_1.shape[0]
    N2 = df_2.shape[0]

    effort_dict = defaultdict(dict)
    for tcr_id in df_2['id']:
        print(tcr_id)
        tcr_row = full_df[full_df['id'] == tcr_id].index.values[0]
        tcr_v_gene = full_df.iloc[tcr_row, :]['v_gene']
        tcr_cdr3 = full_df.iloc[tcr_row, :]['cdr3']
        effort_dict[tcr_id][tcr_v_gene] = {}
        effort_dict[tcr_id][tcr_v_gene][tcr_cdr3] = []

        for trial in range(trial_count):
            df_1_trial, df_2_trial = split_datasets(full_df, N1, N2, tcr_index=tcr_row)

            trial_file_1 = "trial_file_1.csv"
            trial_file_2 = "trial_file_2.csv"

            df_1_trial.iloc[:, [0, 1]].to_csv(trial_file_1, header=False, index=False)
            df_2_trial.iloc[:, [0, 1]].to_csv(trial_file_2, header=False, index=False)

            mass_1, _ = get_mass_objects(df_1_trial, "uniform")
            mass_2, _ = get_mass_objects(df_2_trial, "uniform")

            dist_mat = get_raw_distance_matrix(trial_file_1, trial_file_2, db='/fh/fast/matsen_e/bolson2/transport/iel_data/fake_pubtcrs_db_mouse', exe='../bin/tcrdists', verbose=False)/DMAX
            ot_mat = ot.sinkhorn(mass_1, mass_2, dist_mat, LAMBDA)
            effort_mat = np.multiply(dist_mat, ot_mat)

            efforts = DMAX*N1*effort_mat.sum(axis=0)
            assert efforts.shape[0] == N2


            tcr_position = np.where(df_2_trial['id'] == tcr_id)[0].tolist()[0]

            if method_type is "neighborhood_means":
                dist_mat_df_2 = get_raw_distance_matrix(trial_file_2, trial_file_2, db='/fh/fast/matsen_e/bolson2/transport/iel_data/fake_pubtcrs_db_mouse', exe='../bin/tcrdists', verbose=False)
                ii_nbrhood_mask = (dist_mat_df_2[:, tcr_position] < NEIGHBOR_CUTOFF)
                tcr_effort = np.mean(efforts[ii_nbrhood_mask])
            if method_type is "neighborhood_sums":
                dist_mat_df_2 = get_raw_distance_matrix(trial_file_2, trial_file_2, db='/fh/fast/matsen_e/bolson2/transport/iel_data/fake_pubtcrs_db_mouse', exe='../bin/tcrdists', verbose=False)
                ii_nbrhood_mask = (dist_mat_df_2[:, tcr_position] < NEIGHBOR_CUTOFF)
                tcr_effort = np.sum(efforts[ii_nbrhood_mask])
            elif method_type is "absolute":
                tcr_effort = efforts[tcr_position]
            
            effort_dict[tcr_id][tcr_v_gene][tcr_cdr3].append(tcr_effort)
    return effort_dict


if __name__ == "__main__":
    main_dir = '/home/bolson2/sync/within_gene/'
    file_dir = '/fh/fast/matsen_e/bolson2/transport/iel_data/iels_tcrs_by_mouse/'

    cd4_subject = 'CD4_20'
    dn_subject = 'DN_10' # 797 tcrs
    cd8_subject = 'CD8_15'

    NEIGHBOR_CUTOFF = 48.5
    
    cd4_df = get_processed_df(cd4_subject, file_dir)
    dn_df = get_processed_df(dn_subject, file_dir)
    cd8_df = get_processed_df(cd8_subject, file_dir)
    method_type="neighborhood_means"
    result = do_randomization_test(cd4_df, dn_df, method_type=method_type)
    results_dir = '/home/bolson2/sync/per_tcr/' + method_type
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    with open(results_dir + "/per_tcr.json", 'w') as fp:
        json.dump(result, fp)

