import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.append(os.getcwd())

from common.params import CSV_OUTPUT_DIRNAME, DEFAULT_NEIGHBOR_RADIUS, DIRECTORIES, DIST_MATRICES, JSON_OUTPUT
from python.tcr_clusterer import TCRClusterer

with open(os.path.join(DIRECTORIES[JSON_OUTPUT], 'empirical_fg_bg_nbhd_stats.json')) as f:
    result = json.load(f)

subjects = result.keys()
sample_sizes = {subject: len(result[subject][str(DEFAULT_NEIGHBOR_RADIUS)]) for subject in subjects}
sample_size_threshold = 300
dfs = []
full_dict = {}

for subject in subjects:
    if sample_sizes[subject] > sample_size_threshold:
        subject_distance_matrix = np.loadtxt(os.path.join(DIRECTORIES[DIST_MATRICES], subject + '.csv'), dtype='i', delimiter=',')
        info_dict = result[subject][str(DEFAULT_NEIGHBOR_RADIUS)]
        score_dict = {tcr: tcr_info['foreground']['score'] for tcr, tcr_info in info_dict.items()}

        tcr_clusterer = TCRClusterer(subject_distance_matrix, score_dict)
        
        df = tcr_clusterer.df
        df['subject'] = subject
        dfs.append(df)

        full_dict[subject] = tcr_clusterer.motif_dict

full_df = pd.concat(dfs)

if not os.path.exists(CSV_OUTPUT_DIRNAME):
    os.makedirs(CSV_OUTPUT_DIRNAME)

full_df.to_csv(os.path.join(CSV_OUTPUT_DIRNAME, "motif.csv"), index=False)

with open(os.path.join(DIRECTORIES[JSON_OUTPUT], "motif.json"), "w") as fp:
    json.dump(full_dict, fp)
