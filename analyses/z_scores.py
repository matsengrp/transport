from collections import defaultdict
import json
import os
from random import sample
import sys

import numpy as np
import pandas as pd

sys.path.append('.')

from common.params import DIRECTORIES, DMAX, JSON_OUTPUT
from python.randomization_test import RandomizationTest
from python.tcr_multi_clusterer import TCRMultiClusterer
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

    cd4_filename = get_filename_from_subject(cd4_subject, file_dir)
    dn_filename = get_filename_from_subject(dn_subject, file_dir)
    cd8_filename = get_filename_from_subject(cd8_subject, file_dir)

    obs_scorer = TCRScorer(file_1=cd4_filename, file_2=dn_filename, species="mouse")
    obs_scores = obs_scorer.enrichment_dict

    cd4_df = get_df_from_file(cd4_filename)
    dn_df = get_df_from_file(dn_filename)

    randomization_test = RandomizationTest(cd4_df, dn_df, obs_scores, species="mouse")
    result = randomization_test.enrichment_dict

    results_dir = DIRECTORIES[JSON_OUTPUT]
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    z_scores = {tcr: tcr_info['z_score'] for tcr, tcr_info in result.items()}
    with open(os.path.join(results_dir,  "rand_z_scores.json"), 'w') as fp:
        json.dump(z_scores, fp)

    p_values = {tcr: tcr_info['p_value'] for tcr, tcr_info in result.items()}
    with open(os.path.join(results_dir,  "rand_p_values.json"), 'w') as fp:
        json.dump(z_scores, fp)

    clusterer = TCRMultiClusterer(
        file_1=cd4_filename,
        file_2=dn_filename,
        species="mouse",
        outdir=results_dir
    )
    pd.DataFrame(clusterer.result).transpose().to_csv(os.path.join(results_dir, 'cluster_df.csv'))
