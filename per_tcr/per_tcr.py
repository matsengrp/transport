from random import sample
import sys

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

def do_randomization_test(df_1, df_2, trial_count=10):
    tcr_id = 'DN_11_846' # This TCR has a CDR3 of CALGDH, which is unusally short
    full_df = pd.concat([df_1, df_2], axis=0).reset_index(drop=True)
    tcr_row = full_df[full_df['id'] == tcr_id].index.values[0]
    N1 = df_1.shape[0]
    N2 = df_2.shape[0]

    df_1_trial, df_2_trial = split_datasets(full_df, N1, N2, tcr_index=tcr_row)

if __name__ == "__main__":
    main_dir = '/home/bolson2/sync/within_gene/'
    file_dir = '/fh/fast/matsen_e/bolson2/transport/iel_data/iels_tcrs_by_mouse/'

    cd4_subject = 'CD4_17'
    dn_subject = 'DN_11'
    
    cd4_df = get_processed_df(cd4_subject, file_dir)
    dn_df = get_processed_df(dn_subject, file_dir)

    result = do_randomization_test(cd4_df, dn_df)

