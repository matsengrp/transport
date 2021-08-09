import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.append('.')

from python.tcr_dist import TCRDist
from config import CONFIG

def get_cluster_memberships(query_file, cluster_dict):
    centroid_file = "centroid.tcrs"
    np.savetxt(
        os.path.join(CONFIG["TMP_OUTPUT"], centroid_file), 
        [cluster_dict['centroid']], 
        fmt="%s"
    )
    dist_mat = TCRDist(species="mouse").get_raw_distance_matrix(
        query_file, 
        os.path.join(CONFIG["TMP_OUTPUT"], centroid_file), 
        output_dir=""
    )
    return list(dist_mat[:, 0] < cluster_dict['radius'])

if __name__ == "__main__":
    file_dir = CONFIG["IEL_DATA_DIR"] 
    hmm_output_dir = CONFIG["HMM_CD4_DN_OUTPUT"]
    results_dir = CONFIG["IEL_CLUSTER_OUTPUT"]

    dn_subjects = [f"DN_{x}" for x in range(1, 23 + 1)]
    cd4_subjects = [f"CD4_{x}" for x in range(1, 23 + 1)]
    cd8_subjects = [f"CD8_{x}" for x in range(1, 23 + 1)]
    subjects = dn_subjects + cd4_subjects + cd8_subjects
    for subject in subjects:
        f = os.path.join(file_dir, f"{subject}_B.tcrs")
        cluster_df = pd.read_csv(f, header=None)
        cluster_df['motif'] = "N/A"
        for cluster, motif_name in zip([1, 2, 3], ["Tremont", "Revere", "Ida"]):
            with open(os.path.join(hmm_output_dir,  f"cluster_{cluster}", "cluster_info.json")) as fp:
                cluster_dict = json.load(fp)
            
            memberships = get_cluster_memberships(f, cluster_dict)  
            cluster_df.loc[memberships, 'motif'] = motif_name

        subject_results_dir = os.path.join(results_dir, subject)
        if not os.path.exists(subject_results_dir):
            os.makedirs(subject_results_dir)
        cluster_df.to_csv(os.path.join(subject_results_dir, 'cluster_df.csv'), header=False, index=False)
