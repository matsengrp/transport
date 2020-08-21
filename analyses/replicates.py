## this script compares neighbor-summed per-tcr-efforts (ie row sums of the effort matrix) within a
## set of 'foreground' repertoires and between the foreground repertoires and a set of background 
## repertoires. Here the foreground repertoires are the DN beta repertoires across the mice and the
## background reps are the CD4s (see fg_reptag and bg_reptag below)
##
## The idea is that each foreground tcr will have a distribution of efforts in the comparisons
## to the other foreground repertoires and also a distribution of efforts in the comparisons to
## the background repertoires. For a wonky tcr that's an outlier, both those distributions will
## be large. For a tcr that has more neighbors in the foreground repertoires than in the background
## repertoires, the two distributions will be different and we can attach a significance to that.

from collections import defaultdict
import json

import numpy as np
import ntpath

from glob import glob
import os
from scipy.stats import mannwhitneyu, ttest_ind
import sys

sys.path.append(os.getcwd())
from common.params import DEFAULT_NEIGHBOR_RADIUS, DIRECTORIES
from python.tcr_scorer import TCRScorer

def get_enrichments(scorer, radius):
    scorer.compute_enrichments(neighbor_radius=radius)
    return list(scorer.enrichment_dict.values())

def get_neighbor_counts(scorer):
    return list(scorer.neighbor_counts.values())

# on rhino1:
seq_data_dir = '/loc/no-backup/pbradley/share/pot_data/iels_tcrs_by_mouse/'

for _, d in DIRECTORIES.items():
    if not os.path.exists(d):
        os.makedirs(d)

## fg = foreground, bg = background
fg_reptag = 'CD4'
bg_reptag = 'DN'
chain = 'B'

fg_repfiles = sorted(glob('{}{}_*_{}.tcrs'.format(seq_data_dir, fg_reptag, chain)))
bg_repfiles = sorted(glob('{}{}_*_{}.tcrs'.format(seq_data_dir, bg_reptag, chain)))

neighbor_cutoff_radii = [i + .5 for i in range(0, 200, 10)] # distance at which two single-chain tcrs are considered nbrs

cd4_subject = 'CD4_17_B.tcrs'
dn_subject = 'DN_15_B.tcrs'

file_dir = '/fh/fast/matsen_e/bolson2/transport/iel_data/iels_tcrs_by_mouse/'

cd4_file = os.path.join(file_dir, cd4_subject)


scores_by_cutoff_dict= defaultdict(dict)
z_score_dict = defaultdict(dict)
p_val_dict = defaultdict(dict)
for bg_repfile in bg_repfiles:
    obs_scorer = TCRScorer(file_1=cd4_file, file_2=bg_repfile, species="mouse")
    obs_scores = obs_scorer.enrichment_dict.values()

    tcrs = [x[:-1] for x in open(bg_repfile,'r')]
    unique_tcrs = list(dict.fromkeys(tcrs))
    N1 = len(unique_tcrs)

    fg_scorers = [TCRScorer(x, bg_repfile, species="mouse") for x in fg_repfiles]
    bg_scorers = [TCRScorer(x, bg_repfile, species="mouse") for x in bg_repfiles if x != bg_repfile]

    fg_enrichments = [get_enrichments(fg_scorer, DEFAULT_NEIGHBOR_RADIUS) for fg_scorer in fg_scorers]
    bg_enrichments = [get_enrichments(bg_scorer, DEFAULT_NEIGHBOR_RADIUS) for bg_scorer in bg_scorers]

    T_fg_enrichments = np.array( fg_enrichments ).transpose()
    T_bg_enrichments = np.array( bg_enrichments ).transpose()

    assert T_fg_enrichments.shape == (N1, len(fg_repfiles))
    assert T_bg_enrichments.shape == (N1, len(bg_repfiles) - 1)

    T_fg_means = np.mean(T_fg_enrichments, axis=1)
    T_fg_stddevs = np.std(T_fg_enrichments, axis=1)
    T_bg_means = np.mean(T_bg_enrichments, axis=1)
    T_bg_stddevs = np.std(T_bg_enrichments, axis=1)

    fg_z_scores = [(obs - mean)/std for obs, mean, std in zip(obs_scores, T_fg_means, T_fg_stddevs)]
    bg_z_scores = [(obs - mean)/std for obs, mean, std in zip(obs_scores, T_bg_means, T_bg_stddevs)]

    subject = ntpath.basename(bg_repfile)
    self_distance_matrix = obs_scorer.repertoire_2.distance_matrix
    subject = ntpath.basename(bg_repfile)
    np.savetxt(
        os.path.join(
            DIRECTORIES["dist_matrices"],
            subject + ".csv", 
        ),
        self_distance_matrix,
        delimiter=",",
        fmt='%i'
    )

    z_score_dict[subject] = {tcr: {"foreground": fg_z_score, "background": bg_z_score} for tcr, fg_z_score, bg_z_score in zip(unique_tcrs, fg_z_scores, bg_z_scores)}
    
    result = {"foreground": {"means": T_fg_means.tolist(), "std_devs": T_fg_stddevs.tolist()}, "background": {"means": T_bg_means.tolist(), "std_devs": T_bg_stddevs.tolist()}}
    with open(os.path.join(DIRECTORIES["json_output"], "empirical_fg_bg_stats.json"), "w") as fp:
        json.dump(result, fp)

    for ii in range(N1):

        _, mann_whitney_pval = mannwhitneyu(T_bg_enrichments[ii, :], T_fg_enrichments[ii, :])
        _, t_test_p_value = ttest_ind(T_bg_enrichments[ii, :], T_fg_enrichments[ii, :])

        # crude multiple-testing correction:
        mann_whitney_pval *= N1
        t_test_p_value *= N1

        p_val_dict[subject][unique_tcrs[ii]] = mann_whitney_pval

    nbhd_means_by_cutoff = defaultdict(dict)
    cutoff_means = defaultdict(dict)
    cutoff_sds = defaultdict(dict)
    for neighbor_cutoff_radius in neighbor_cutoff_radii:
        bg_enrichments = np.array([get_enrichments(bg_scorer, radius=neighbor_cutoff_radius) for bg_scorer in bg_scorers]).transpose()
        fg_enrichments = np.array([get_enrichments(fg_scorer, radius=neighbor_cutoff_radius) for fg_scorer in fg_scorers]).transpose()
        bg_neighbor_counts = np.array([get_neighbor_counts(bg_scorer) for bg_scorer in bg_scorers]).transpose()
        fg_neighbor_counts = np.array([get_neighbor_counts(fg_scorer) for fg_scorer in fg_scorers]).transpose()

        cutoff_means[neighbor_cutoff_radius] = defaultdict(dict)
        cutoff_sds[neighbor_cutoff_radius] = defaultdict(dict)
        for ii in range(N1):
            cutoff_means[neighbor_cutoff_radius][unique_tcrs[ii]] = {
                    "foreground": {"score": np.mean(fg_enrichments[ii, :]), "neighbor_count": np.mean(fg_neighbor_counts[ii, :])},
                    "background": {"score": np.mean(bg_enrichments[ii, :]), "neighbor_count": np.mean(bg_neighbor_counts[ii, :])}
            }

            cutoff_sds[neighbor_cutoff_radius][unique_tcrs[ii]] = {
                "foreground": np.std(fg_enrichments[ii, :]),
                "background": np.std(bg_enrichments[ii, :])
            }

    scores_by_cutoff_dict[subject] = cutoff_means

with open(os.path.join(DIRECTORIES["json_output"], "empirical_fg_bg_nbhd_stats.json"), "w") as fp:
    json.dump(scores_by_cutoff_dict, fp)

with open(os.path.join(DIRECTORIES["json_output"], "replicate_z_scores.json"), "w") as fp:
    json.dump(z_score_dict, fp)
