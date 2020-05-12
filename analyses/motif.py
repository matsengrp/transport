from collections import defaultdict
from itertools import compress
import json
import operator
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize

sys.path.append(os.getcwd())

from common.params import CSV_OUTPUT_DIRNAME, DEFAULT_NEIGHBOR_RADIUS, DIRECTORIES, DIST_MATRICES, JSON_OUTPUT, TMP_OUTPUT

with open(os.path.join(DIRECTORIES[JSON_OUTPUT], 'empirical_fg_bg_nbhd_stats.json')) as f:
    result = json.load(f)

subjects = result.keys()
sample_sizes = {subject: len(result[subject][str(DEFAULT_NEIGHBOR_RADIUS)]) for subject in subjects}
sample_size_threshold = 300
dfs = []

full_dict = {}
for subject in subjects:
    if sample_sizes[subject] > sample_size_threshold:
        subject_distance_matrix = np.loadtxt(os.path.join(DIRECTORIES[DIST_MATRICES], subject + '.csv'), dtype='i', delimiter=',')
        info_dict = result[subject][str(DEFAULT_NEIGHBOR_RADIUS)]
        score_dict = {tcr: tcr_info['foreground']['score'] for tcr, tcr_info in info_dict.items()}
        unique_tcrs = list(dict.fromkeys(score_dict.keys()))
        scores = list(score_dict.values())
        
        max_score_tcr = max(score_dict.items(), key=operator.itemgetter(1))[0]
        index = unique_tcrs.index(max_score_tcr)
        
        enrichment_threshold = np.quantile(scores, 0.75)
        enrichment_mask = (scores > enrichment_threshold)
        
        alignment_infile = os.path.join(DIRECTORIES[TMP_OUTPUT], "cluster_cdr3s.fasta")
        alignment_outfile = os.path.join(DIRECTORIES[TMP_OUTPUT], "cluster_cdr3s.sto")
        hmm_filename = "TRB_mouse.hmm"
        motif_hmm_file = os.path.join(DIRECTORIES[TMP_OUTPUT], "motif.hmm")
        motif_hmm_stats_file = os.path.join(DIRECTORIES[TMP_OUTPUT], "motif_stats.txt")
        
        if not os.path.exists(hmm_filename):
            os.mknod(hmm_filename)
        
        radii = [i + .5 for i in range(0, 200, 1)]
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
        
            records = [SeqRecord(Seq(cdr3), id=tcr) for tcr, cdr3 in zip(neighborhood_tcrs, neighborhood_cdr3s)]
            with open(alignment_infile, "w") as output_handle:
                SeqIO.write(records, output_handle, "fasta")
            
            os.system('hmmalign {} {} > {}'.format(hmm_filename, alignment_infile, alignment_outfile))
            os.system('hmmbuild {} {}'.format(motif_hmm_file, alignment_outfile))
            os.system('hmmstat {} > {}'.format(motif_hmm_file, motif_hmm_stats_file))

            fields = ['idx', 'name', 'accession', 'nseq', 'eff_nseq', 'M', 'relent', 'info', 'p relE', 'compKL']
            with open(motif_hmm_stats_file, 'r') as f:
                for line in f:
                    print(line)
                    if line.startswith('1'):
                        values = line.split()
                        hmm_stats = {field: value for field, value in zip(fields, values)}
            motif_dict[radius]['entropy'] = hmm_stats['relent']
        
        df = pd.DataFrame.from_dict(
            {
                radius: {
                    'cluster_size': v['cluster_size'],
                    'mean_enrichment': v['mean_enrichment'],
                    'entropy': v['entropy'],
                } for radius, v in motif_dict.items()
            }
        )
        df = df.transpose()
        df['radius'] = motif_dict.keys()
        df['subject'] = subject
        dfs.append(df)

    full_dict[subject] = motif_dict

full_df = pd.concat(dfs)

if not os.path.exists(CSV_OUTPUT_DIRNAME):
    os.makedirs(CSV_OUTPUT_DIRNAME)

full_df.to_csv(os.path.join(CSV_OUTPUT_DIRNAME, "motif.csv"), index=False)

with open(os.path.join(DIRECTORIES[JSON_OUTPUT], "motif.json"), "w") as fp:
    json.dump(full_dict, fp)
