import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.append('.')

from common.params import DIRECTORIES, TMP_OUTPUT
from python.tcr_dist import TCRDist

def get_cluster_memberships(query_file, cluster_dict):
    centroid_file = "centroid.tcrs"
    np.savetxt(
        os.path.join(DIRECTORIES[TMP_OUTPUT], centroid_file), 
        [cluster_dict['centroid']], 
        fmt="%s"
    )
    dist_mat = TCRDist(species="mouse").get_raw_distance_matrix(
        query_file, 
        os.path.join(DIRECTORIES[TMP_OUTPUT], centroid_file)
    )
    return list(dist_mat[:, 0] < cluster_dict['radius'])

if __name__ == "__main__":
    file_dir = '/fh/fast/matsen_e/bolson2/transport/iel_data/iels_tcrs_by_mouse/'

    hmm_output_dir = "output/hmm/cd4_dn"

    results_dir = "output/iel_clusters"
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    subjects = [f"DN_{x}" for x in range(1, 23 + 1)]
    for subject in subjects:
        f = os.path.join(file_dir, f"{subject}_B.tcrs")
        cluster_df = pd.read_csv(f, header=None)
        cluster_df['motif'] = "N/A"
        for cluster, motif_name in zip([1, 2, 3], ["Ida", "Revere", "Tremont"]):
            with open(os.path.join(hmm_output_dir,  f"cluster_{cluster}", "cluster_info.json")) as fp:
                cluster_dict = json.load(fp)
            
            memberships = get_cluster_memberships(f, cluster_dict)  
            cluster_df.loc[memberships, 'motif'] = motif_name

        subject_results_dir = os.path.join(results_dir, subject)
        if not os.path.exists(subject_results_dir):
            os.makedirs(subject_results_dir)
        cluster_df.to_csv(os.path.join(subject_results_dir, 'cluster_df.csv'))
