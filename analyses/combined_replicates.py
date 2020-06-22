from glob import glob
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

sys.path.append(os.getcwd())
from common.params import CSV_OUTPUT_DIRNAME, DIRECTORIES, IEL_DATA_DIR, JSON_OUTPUT, TMP_OUTPUT
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

scorer = TCRScorer(file_1=cd4_file, file_2=dn_file)
dn_self_dist_mat = scorer.repertoire_2.distance_matrix
tcr_clusterer = TCRClusterer(self_distance_matrix=dn_self_dist_mat, score_dict=scorer.enrichment_dict)
tcr_clusterer.cluster_df.to_csv(os.path.join(CSV_OUTPUT_DIRNAME, "combined_cluster_df.csv"))

np.savetxt(
    os.path.join(
        DIRECTORIES["dist_matrices"],
        "full_dn.csv",
    ),
    dn_self_dist_mat,
    delimiter=",",
    fmt="%i"
)
