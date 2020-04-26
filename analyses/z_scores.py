from collections import defaultdict
import json
import os
from random import sample
import sys

import numpy as np

sys.path.append('.')

from python.randomization import do_randomization_test, get_filename_from_subject, get_processed_df, split_datasets
from python.utils import get_effort_scores

LAMBDA = 0.01
DMAX = 200

if __name__ == "__main__":
    main_dir = '/home/bolson2/sync/within_gene/'
    file_dir = '/fh/fast/matsen_e/bolson2/transport/iel_data/iels_tcrs_by_mouse/'

    cd4_subject = 'CD4_16'
    dn_subject = 'DN_18'
    cd8_subject = 'CD8_15'

    NEIGHBOR_CUTOFF = 50.5
    
    cd4_df = get_processed_df(cd4_subject, file_dir)
    dn_df = get_processed_df(dn_subject, file_dir)
    cd8_df = get_processed_df(cd8_subject, file_dir)
    method_type="neighborhood_sums"

    obs_scores = get_effort_scores(cd4_df, dn_df)

    result = do_randomization_test(cd4_df, dn_df, method_type=method_type)

    # Downsample all score distributions to smallest observed count
    min_sample_size = np.min([len(x) for x in result.values()])
    for tcr, scores in result.items():
        result[tcr] = sample(scores, min_sample_size)

    # Finally, compute the mean score for each TCR:
    mean_result = {tcr: np.mean(scores) for tcr, scores in result.items()}

    results_dir = '/home/bolson2/sync/per_tcr/' + method_type
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    with open(results_dir + "/per_tcr.json", 'w') as fp:
        json.dump(mean_result, fp)


    z_scores = defaultdict()
    for obs_score, sim_scores, tcr in zip(obs_scores, result.values(), dn_df['TCR']):
        z_scores[tcr] = (obs_score - np.mean(sim_scores))/np.std(sim_scores)
    with open(results_dir + "/z_scores.json", 'w') as fp:
        json.dump(z_scores, fp)
