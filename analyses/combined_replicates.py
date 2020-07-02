from glob import glob
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.append(os.getcwd())
from common.params import CSV_OUTPUT_DIRNAME, DIRECTORIES, HMM_OUTPUT, IEL_DATA_DIR, JSON_OUTPUT, TMP_OUTPUT
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

def run_clustering_step(file_1, file_2, cluster, score_dict=None):
    scorer = TCRScorer(file_1=file_1, file_2=file_2)
    if score_dict is None:
        score_dict=scorer.enrichment_dict
    rep_2_self_dist_mat = scorer.repertoire_2.distance_matrix
    tcr_clusterer = TCRClusterer(self_distance_matrix=rep_2_self_dist_mat, score_dict=score_dict, cluster_label=f"cluster_{cluster}")
    return scorer, tcr_clusterer

all_clusters = []

initial_scorer = TCRScorer(file_1=cd4_file, file_2=dn_file)
initial_rep_2_self_dist_mat = initial_scorer.repertoire_2.distance_matrix
initial_score_dict = initial_scorer.enrichment_dict
initial_clusterer = TCRClusterer(self_distance_matrix=initial_rep_2_self_dist_mat, score_dict=initial_score_dict, cluster_label="cluster_1")
result = {tcr: {'score': score, 'cluster': 0} for tcr, score in initial_scorer.enrichment_dict.items()} 
for tcr in initial_clusterer.cluster_dict['tcrs']:
    result[tcr]['cluster'] = 1

sub_repertoire_tcrs = [tcr for tcr in initial_scorer.repertoire_2.unique_tcrs if tcr not in initial_clusterer.cluster_dict['tcrs']]
hmmer_manager = HMMerManager()
if not os.path.exists(DIRECTORIES[HMM_OUTPUT]):
    os.makedirs(DIRECTORIES[HMM_OUTPUT])

hmmer_manager.build_hmm_from_sequences(
    [s.split(',')[1] for s in initial_clusterer.cluster_dict['tcrs']],
    hmm_filename=os.path.join(DIRECTORIES[HMM_OUTPUT], 'cluster_1.hmm'),
    alignment_outfilename=os.path.join(DIRECTORIES[HMM_OUTPUT], 'cluster_1.sto'),
    fasta_filename=os.path.join(DIRECTORIES[HMM_OUTPUT], 'cluster_1.fasta'),
    plot_filename=os.path.join(DIRECTORIES[HMM_OUTPUT], 'cluster_1.png')
)
cluster = 2
sub_dn_file = os.path.join(DIRECTORIES[TMP_OUTPUT], 'sub_dn.csv')
while cluster < 10:
    np.savetxt(sub_dn_file, sub_repertoire_tcrs, fmt="%s")
    score_dict = {k: initial_score_dict[k] for k in sub_repertoire_tcrs}
    current_scorer, current_clusterer = run_clustering_step(file_1=cd4_file, file_2=sub_dn_file, cluster=cluster, score_dict=score_dict) 
    sub_repertoire_tcrs = [tcr for tcr in current_scorer.repertoire_2.unique_tcrs if tcr not in current_clusterer.cluster_dict['tcrs']]
    for tcr in current_clusterer.cluster_dict['tcrs']:
        result[tcr]['cluster'] = cluster
    hmmer_manager = HMMerManager()
    hmmer_manager.build_hmm_from_sequences(
        [s.split(',')[1] for s in current_clusterer.cluster_dict['tcrs']],
        hmm_filename=os.path.join(DIRECTORIES[HMM_OUTPUT], 'cluster_{}.hmm'.format(cluster)),
        alignment_outfilename=os.path.join(DIRECTORIES[HMM_OUTPUT], 'cluster_{}.sto'.format(cluster)),
        fasta_filename=os.path.join(DIRECTORIES[HMM_OUTPUT], 'cluster_{}.fasta'.format(cluster)),
        plot_filename=os.path.join(DIRECTORIES[HMM_OUTPUT], 'cluster_{}.png'.format(cluster))
    )
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


