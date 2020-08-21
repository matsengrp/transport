from collections import defaultdict
import json
import os
from random import sample
import sys

import numpy as np
import pandas as pd

sys.path.append('.')

from common.params import DIRECTORIES, DMAX, JSON_OUTPUT, TMP_OUTPUT
from python.hmmer_manager import HMMerManager
from python.randomization_test import RandomizationTest
from python.tcr_multi_clusterer import TCRMultiClusterer
from python.tcr_scorer import TCRScorer
from python.utils import extract_cdr3s, get_df_from_file

def get_filename_from_subject(subject, file_dir):
    filename = file_dir + subject + '_B.tcrs'
    return filename


if __name__ == "__main__":
    file_dir = '/fh/fast/matsen_e/bolson2/transport/iel_data/iels_tcrs_by_mouse/'

    cd4_subject = 'CD4_17'
    dn_subject = 'DN_15'
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
    cluster_df = pd.DataFrame(clusterer.result).transpose()
    cluster_df['tcr'] = cluster_df.index
    cluster_df['motif'] = "N/A"

    hmm_output_dir = "output/hmm/cd4_dn"
    hmmer_manager = HMMerManager()
    sequences = list(clusterer.result.keys())
    cdr3s = extract_cdr3s(sequences)

    for cluster, motif_name in zip([1, 2, 3], ["Ida", "Revere", "Tremont"]):
        search_result = hmmer_manager.run_hmmsearch(
            hmm_filename=os.path.join(hmm_output_dir, f"cluster_{cluster}/seqs.hmm"),
            query_sequences=extract_cdr3s(list(clusterer.result.keys())),
            outfile=os.path.join(DIRECTORIES[TMP_OUTPUT], f"{motif_name}_hmmsearch.out"),
            sequence_ids=sequences
        )
        e_values = [float(i['e_value']) for i in search_result]
        motifs = cluster_df['motif']
        for i in range(len(e_values)):
            if e_values[i] < 1e-8 and motifs[i] == "N/A":
                tcr = search_result[i]['target_name']
                cluster_df.loc[cluster_df['tcr'] == tcr, 'motif'] = motif_name

    cluster_df.to_csv(os.path.join(results_dir, 'cluster_df.csv'))
    import pdb; pdb.set_trace()
