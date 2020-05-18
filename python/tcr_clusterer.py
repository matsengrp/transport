from collections import defaultdict
from itertools import compress
import operator
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pandas as pd

from common.params import DEFAULT_NEIGHBOR_RADIUS, DIRECTORIES, TMP_OUTPUT

class TCRClusterer():
    def __init__(self, self_distance_matrix, score_dict, radius=DEFAULT_NEIGHBOR_RADIUS):
        unique_tcrs = list(dict.fromkeys(score_dict.keys()))
        scores = list(score_dict.values())

        max_score_tcr = max(score_dict.items(), key=operator.itemgetter(1))[0]
        index = unique_tcrs.index(max_score_tcr)

        enrichment_threshold = np.quantile(scores, 0.5)
        enrichment_mask = (scores > enrichment_threshold)

        alignment_infile = os.path.join(DIRECTORIES[TMP_OUTPUT], "cluster_cdr3s.fasta")
        alignment_outfile = os.path.join(DIRECTORIES[TMP_OUTPUT], "cluster_cdr3s.sto")
        hmm_filename = "TRB_mouse.hmm"
        motif_hmm_file = os.path.join(DIRECTORIES[TMP_OUTPUT], "motif.hmm")
        motif_hmm_stats_file = os.path.join(DIRECTORIES[TMP_OUTPUT], "motif_stats.txt")
        
        if not os.path.exists(hmm_filename):
            os.mknod(hmm_filename)
        
        radii = [i + .5 for i in range(0, 200, 5)]
        self.motif_dict = defaultdict(dict)
        previous_radius = -1
        for radius in radii:
            neighborhood_mask = (self_distance_matrix[index, :] < radius)
            annulus_mask = (self_distance_matrix[index, :] < radius) & (self_distance_matrix[index, :] > previous_radius) & enrichment_mask
            full_mask = neighborhood_mask & enrichment_mask
            neighborhood_enrichments = list(compress(scores, full_mask))
            neighborhood_tcrs = list(compress(unique_tcrs, full_mask))
            neighborhood_cdr3s = [s.split(',')[1] for s in neighborhood_tcrs]
            annulus_enrichments = list(compress(scores, annulus_mask))
            self.motif_dict[radius] = {
                "mean_enrichment": np.mean(neighborhood_enrichments),
                "annulus_enrichment": np.mean(annulus_enrichments) if annulus_enrichments != [] else None,
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
            self.motif_dict[radius]['entropy'] = hmm_stats['relent']

            previous_radius = radius
        
        self.df = pd.DataFrame.from_dict(
            {
                radius: {
                    'cluster_size': v['cluster_size'],
                    'mean_enrichment': v['mean_enrichment'],
                    'annulus_enrichment': v['annulus_enrichment'],
                    'entropy': v['entropy'],
                } for radius, v in self.motif_dict.items()
            }
        )
        self.df = self.df.transpose()
        self.df['radius'] = self.motif_dict.keys()
