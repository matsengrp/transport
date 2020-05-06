from collections import defaultdict
from itertools import compress
import json
import operator
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize

sys.path.append(os.getcwd())

from common.params import CSV_OUTPUT_DIRNAME, DEFAULT_NEIGHBOR_RADIUS, DIRECTORIES, DIST_MATRICES, JSON_OUTPUT, TMP_OUTPUT

with open(os.path.join(DIRECTORIES[JSON_OUTPUT], 'empirical_fg_bg_nbhd_stats.json')) as f:
    result = json.load(f)

subject = 'DN_10_B.tcrs'
subject_distance_matrix = np.loadtxt(os.path.join(DIRECTORIES[DIST_MATRICES], subject + '.csv'), dtype='i', delimiter=',')
info_dict = result[subject][str(DEFAULT_NEIGHBOR_RADIUS)]
score_dict = {tcr: tcr_info['foreground']['score'] for tcr, tcr_info in info_dict.items()}
unique_tcrs = list(dict.fromkeys(score_dict.keys()))
scores = list(score_dict.values())

max_score_tcr = max(score_dict.items(), key=operator.itemgetter(1))[0]
index = unique_tcrs.index(max_score_tcr)

enrichment_threshold = np.quantile(scores, 0.5)
enrichment_mask = (scores > enrichment_threshold)

txt_filename = os.path.join(DIRECTORIES[TMP_OUTPUT], "cluster_cdr3s.txt")
hmm_filename = os.path.join(DIRECTORIES[TMP_OUTPUT], "alignment.hmm")

if not os.path.exists(hmm_filename):
    os.mknod(hmm_filename)

radii = [i + .5 for i in range(0, 100, 1)]
motif_dict = defaultdict(dict)
for radius in radii:
    neighborhood_mask = (subject_distance_matrix[index, :] < radius)
    full_mask = neighborhood_mask & enrichment_mask
    neighborhood_enrichments = list(compress(scores, full_mask))
    neighborhood_tcrs = list(compress(unique_tcrs, full_mask))
    neighborhood_cdr3s = [s.split(',')[1] for s in neighborhood_tcrs]
    motif_dict[radius] = {
        "mean_enrichment": np.mean(neighborhood_enrichments),
        "tcrs": neighborhood_tcrs,
        "cluster_size": int(full_mask.sum())
    }
    print("{}, {}".format(neighborhood_mask.sum(), np.mean(neighborhood_enrichments)))

    np.savetxt(
        txt_filename,
        neighborhood_cdr3s,
        fmt="%s"
    )
    #os.system('hmmbuild {} {}'.format(hmm_filename, txt_filename))

df = pd.DataFrame.from_dict({radius: {'cluster_size': v['cluster_size'], 'mean_enrichment': v['mean_enrichment']} for radius, v in motif_dict.items()})
df = df.transpose()
df['radius'] = motif_dict.keys()

if not os.path.exists(CSV_OUTPUT_DIRNAME):
    os.makedirs(CSV_OUTPUT_DIRNAME)

df.to_csv(os.path.join(CSV_OUTPUT_DIRNAME, "motif.csv"), index=False)

with open(os.path.join(DIRECTORIES[JSON_OUTPUT], "motif.json"), "w") as fp:
    json.dump(motif_dict, fp)
