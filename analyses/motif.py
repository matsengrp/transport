from itertools import compress
import json
import operator
import os
import sys

import numpy as np

sys.path.append(os.getcwd())

from common.params import DEFAULT_NEIGHBOR_RADIUS, DIRECTORIES, DIST_MATRICES, JSON_OUTPUT

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

radii = [i + .5 for i in range(0, 100, 5)]
for radius in radii:
    neighborhood_mask = (subject_distance_matrix[index, :] < radius)
    neighborhood_enrichments = list(compress(scores, neighborhood_mask & enrichment_mask))
    print("{}, {}".format(neighborhood_mask.sum(), np.mean(neighborhood_enrichments)))
