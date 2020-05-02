from collections import Counter
import os
from random import sample

import numpy as np
import pandas as pd

from python.tcr_scorer import TCRScorer
from python.utils import append_id_column, get_df_from_file

def get_filename_from_subject(subject, file_dir):
    filename = file_dir + subject + '_B.tcrs'
    return filename

def get_processed_df(subject, file_dir):
    filename = get_filename_from_subject(subject, file_dir)
    df = get_df_from_file(filename)
    df = append_id_column(df)
    return df

def split_datasets(full_df, N1, N2, filename_1, filename_2):
    indices_to_split = full_df.index.tolist()
    df_1_trial_indices = sample(indices_to_split, N1)
    df_1_trial = full_df.iloc[df_1_trial_indices]
    df_2_trial = full_df.iloc[~full_df.index.isin(df_1_trial_indices)]

    df_1_trial.iloc[:, [0, 1]].to_csv(filename_1, header=False, index=False)
    df_2_trial.iloc[:, [0, 1]].to_csv(filename_2, header=False, index=False)

    return df_1_trial, df_2_trial

def do_randomization_test(df_1, df_2, trial_count=100, output_dir="tmp_output"):
    full_df = pd.concat([df_1, df_2], axis=0).reset_index(drop=True)
    N1 = df_1.shape[0]
    N2 = df_2.shape[0]
    df_2_tcrs = list(dict.fromkeys(df_2['tcr']))
    effort_dict = {tcr: [] for tcr in df_2_tcrs}

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for trial in range(trial_count):
        trial_1_file = os.path.join(output_dir, "trial_1.csv")
        trial_2_file = os.path.join(output_dir, "trial_2.csv")

        df_1_trial, df_2_trial = split_datasets(full_df, N1, N2, trial_1_file, trial_2_file)

        trial_scorer = TCRScorer(file_1=trial_1_file, file_2=trial_2_file)

        for tcr, score in trial_scorer.effort_dict.items():
            if tcr in df_2_tcrs:
                effort_dict[tcr].append(score)

    return effort_dict
