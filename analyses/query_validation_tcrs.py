import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.append(os.getcwd())
from python.tcr_dist import TCRDist

tcr_dir = "data/R"
yfv_tcr_file = os.path.join(tcr_dir, "yfv_tcrs.csv")
cmv_tcr_file = os.path.join(tcr_dir, "cmv_tcrs.csv")

clusters = ["cluster_" + str(i) for i in range(1, 11)]
hit_df = pd.DataFrame(columns=['antigen', 'comparison', 'cluster', 'hits'])

comparisons = [
        "P1_0_-7", "P1_0_15", "P1_0_45",
        "P2_0_-7", "P2_0_15", "P2_0_45",

        "Q1_0_-7", "Q1_0_15", "Q1_0_45",
        "Q2_0_-7", "Q2_0_15", "Q2_0_45",

        "S1_0_-7", "S1_0_15", "S1_0_45",
        "S2_0_-7", "S2_0_15", "S2_0_45",
]
for comparison in comparisons:
    for cluster in clusters:
        comparison_dir = os.path.join("output/hmm", comparison)
        cluster_dir = os.path.join(comparison_dir, cluster)
        cluster_info_file = os.path.join(cluster_dir, "cluster_info.json")
        
        with open(cluster_info_file) as f:
            cluster_info = json.load(f)
        centroid = cluster_info['centroid']
        radius = cluster_info['radius']
        
        centroid_file = os.path.join(cluster_dir, "centroid.txt")
        np.savetxt(centroid_file, [centroid], fmt="%s")
        
        yfv_dist_mat = TCRDist(species="human").get_raw_distance_matrix(yfv_tcr_file, centroid_file, output_dir="")
        yfv_hit_count = (yfv_dist_mat < radius).sum()
        cmv_dist_mat = TCRDist(species="human").get_raw_distance_matrix(cmv_tcr_file, centroid_file, output_dir="")
        cmv_hit_count = (cmv_dist_mat < radius).sum()
        hit_df = hit_df.append({'antigen': 'cmv', 'comparison': comparison, 'cluster': cluster, 'hits': cmv_hit_count}, ignore_index=True)
        hit_df = hit_df.append({'antigen': 'yfv', 'comparison': comparison, 'cluster': cluster, 'hits': yfv_hit_count}, ignore_index=True)

hit_df.to_csv('output/hmm/yfv_hits.csv', index=False)
