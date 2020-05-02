from collections import defaultdict
import json
import os
from random import sample
import sys

import numpy as np

sys.path.append('.')

from common.params import DIRECTORIES, DMAX, JSON_OUTPUT
from python.randomization_test import RandomizationTest
from python.tcr_scorer import TCRScorer
from python.utils import get_df_from_file

def get_filename_from_subject(subject, file_dir):
    filename = file_dir + subject + '_B.tcrs'
    return filename

if __name__ == "__main__":
    file_dir = '/fh/fast/matsen_e/bolson2/transport/iel_data/iels_tcrs_by_mouse/'

    cd4_subject = 'CD4_16'
    dn_subject = 'DN_18'
    cd8_subject = 'CD8_15'

    NEIGHBOR_CUTOFF = 50.5
    
    cd4_filename = get_filename_from_subject(cd4_subject, file_dir)
    dn_filename = get_filename_from_subject(dn_subject, file_dir)
    cd8_filename = get_filename_from_subject(cd8_subject, file_dir)

    obs_scorer = TCRScorer(file_1=cd4_filename, file_2=dn_filename)
    obs_scores = obs_scorer.effort_dict

    cd4_df = get_df_from_file(cd4_filename)
    dn_df = get_df_from_file(dn_filename)

    randomization_test = RandomizationTest(cd4_df, dn_df)
    result = randomization_test.effort_dict

    # Downsample all score distributions to smallest observed count
    min_sample_size = np.min([len(x) for x in result.values()])
    for tcr, scores in result.items():
        result[tcr] = sample(scores, min_sample_size)

    # Finally, compute the mean score for each TCR:
    mean_result = {tcr: np.mean(scores) for tcr, scores in result.items()}

    results_dir = DIRECTORIES[JSON_OUTPUT]
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    with open(os.path.join(results_dir, "per_tcr.json"), 'w') as fp:
        json.dump(mean_result, fp)


    z_scores = defaultdict()
    dn_tcrs = list(dict.fromkeys(dn_df['tcr']))
    for obs_score, sim_scores, tcr in zip(obs_scores.values(), result.values(), dn_tcrs):
        z_scores[tcr] = (obs_score - np.mean(sim_scores))/np.std(sim_scores)
    with open(os.path.join(results_dir,  "rand_z_scores.json"), 'w') as fp:
        json.dump(z_scores, fp)
