from collections import defaultdict
import json

import numpy as np
import ntpath

from glob import glob
import os
from scipy.stats import mannwhitneyu, ttest_ind
import sys

sys.path.append(os.getcwd())
from python.tcr_clusterer import TCRClusterer
from python.tcr_scorer import TCRScorer

seq_data_dir = 'data/yfv'

file_1 = os.path.join(seq_data_dir, "P1_0_F1_.txt.top1000.tcrs")
file_2 = os.path.join(seq_data_dir, "P1_15_F1_.txt.top1000.tcrs")
scorer = TCRScorer(file_1=file_1, file_2=file_2, species="human")
rep_2_self_dist_mat = scorer.repertoire_2.distance_matrix
rep_2_self_dist_mat = TCRClusterer(self_distance_matrix=rep_2_self_dist_mat, score_dict=scorer.enrichment_dict)
