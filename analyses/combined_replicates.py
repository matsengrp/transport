from glob import glob
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.append(os.getcwd())
from common.params import CSV_OUTPUT_DIRNAME, DIRECTORIES, IEL_DATA_DIR, JSON_OUTPUT, TMP_OUTPUT
from python.hmmer_manager import HMMerManager
from python.tcr_clusterer import TCRClusterer
from python.tcr_dist import TCRDist
from python.tcr_scorer import TCRScorer
from python.utils import get_df_from_file

def write_full_replicate_dataset(output_filename, reptag, chain="B"):
    repfiles = sorted(glob('{}{}_*_{}.tcrs'.format(IEL_DATA_DIR, reptag, chain)))
    dfs = []
    for repfile in repfiles:
        dfs.append(get_df_from_file(repfile))
    full_df = pd.concat(dfs)
    full_df.iloc[:, [0, 1]].to_csv(output_filename, index=False, header=False)

cd4_reptag = "CD4"
dn_reptag = "DN"
chain = "B"

dn_filename = "all_dn.csv"
cd4_filename = "all_cd4.csv"
dn_file = os.path.join(DIRECTORIES[TMP_OUTPUT], dn_filename)
cd4_file = os.path.join(DIRECTORIES[TMP_OUTPUT], cd4_filename)
write_full_replicate_dataset(dn_file, "DN")
write_full_replicate_dataset(cd4_file, "CD4")

def run_clustering_step(file_1, file_2):
    scorer = TCRScorer(file_1=file_1, file_2=file_2)
    rep_2_self_dist_mat = scorer.repertoire_2.distance_matrix
    tcr_clusterer = TCRClusterer(self_distance_matrix=rep_2_self_dist_mat, score_dict=scorer.enrichment_dict)
    return scorer, tcr_clusterer

initial_scorer = TCRScorer(file_1=cd4_file, file_2=dn_file)
initial_rep_2_self_dist_mat = initial_scorer.repertoire_2.distance_matrix
initial_tcr_clusterer = TCRClusterer(self_distance_matrix=initial_rep_2_self_dist_mat, score_dict=initial_scorer.enrichment_dict)
result = {tcr: {'score': score, 'cluster': 0} for tcr, score in initial_scorer.enrichment_dict.items()} 
for tcr in initial_tcr_clusterer.cluster_dict['tcrs']:
    result[tcr]['cluster'] = 1

sub_repertoire_tcrs = [tcr for tcr in initial_scorer.repertoire_2.unique_tcrs if tcr not in initial_tcr_clusterer.cluster_dict['tcrs']]
cluster = 2
sub_cd4_file = os.path.join(DIRECTORIES[TMP_OUTPUT], 'sub_cd4.csv')
while cluster < 6:
    np.savetxt(sub_cd4_file, sub_repertoire_tcrs, fmt="%s")
    current_scorer, current_clusterer = run_clustering_step(file_1=cd4_file, file_2=sub_cd4_file) 
    sub_repertoire_tcrs = [tcr for tcr in current_scorer.repertoire_2.unique_tcrs if tcr not in current_clusterer.cluster_dict['tcrs']]
    for tcr in current_clusterer.cluster_dict['tcrs']:
        result[tcr]['cluster'] = cluster
    cluster = cluster + 1


pd.DataFrame(result).transpose().to_csv(os.path.join(DIRECTORIES[TMP_OUTPUT], 'multiple_clusters.csv'))

#np.savetxt(
#    os.path.join(
#        DIRECTORIES["dist_matrices"],
#        "full_dn.csv",
#    ),
#    dn_self_dist_mat,
#    delimiter=",",
#    fmt="%i"
#)

sub_repertoire_tcrs = [tcr for tcr in scorer.repertoire_2.unique_tcrs if tcr not in tcr_clusterer.cluster_dict['tcrs']]
sub_dn_file = os.path.join(DIRECTORIES[TMP_OUTPUT], "sub_dn.csv")
scorer = TCRScorer(file_1=cd4_file, file_2=sub_dn_file)


#hmmer_manager = HMMerManager()
#hmmer_manager.build_hmm_from_sequences([s.split(',')[1] for s in tcr_clusterer.cluster_dict['tcrs']])
