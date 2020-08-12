from collections import defaultdict
from itertools import compress
import operator
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pandas as pd

from common.params import CSV_OUTPUT_DIRNAME, DEFAULT_NEIGHBOR_RADIUS, DIRECTORIES, TMP_OUTPUT
from python.hmmer_manager import HMMerManager

class TCRClusterer():
    seg_csv_file = os.path.join(CSV_OUTPUT_DIRNAME, "seg.csv")

    def __init__(self, self_distance_matrix, score_dict, radius=DEFAULT_NEIGHBOR_RADIUS, seg_csv_outfile=None, cluster_label='tmp', outdir=DIRECTORIES[TMP_OUTPUT]):
        if seg_csv_outfile is None:
            seg_csv_outfile = self.seg_csv_file

        self.unique_tcrs = list(dict.fromkeys(score_dict.keys()))
        self.scores = list(score_dict.values())

        max_score_tcr = max(score_dict.items(), key=operator.itemgetter(1))[0]
        max_score_tcr_index = self.unique_tcrs.index(max_score_tcr)

        enrichment_threshold = np.quantile(self.scores, 0.5)
        enrichment_mask = (self.scores > enrichment_threshold)

        
        radii = [i + .5 for i in range(0, 200, 5)]
        self.all_radii_dict = defaultdict(dict)
        previous_radius = -1
        hmmer_manager = HMMerManager()

        for radius in radii:
            neighborhood_mask = (self_distance_matrix[max_score_tcr_index, :] < radius)
            annulus_mask = (self_distance_matrix[max_score_tcr_index, :] < radius) & (self_distance_matrix[max_score_tcr_index, :] > previous_radius) & enrichment_mask
            full_mask = neighborhood_mask & enrichment_mask
            neighborhood_enrichments = list(compress(self.scores, full_mask))
            neighborhood_tcrs = list(compress(self.unique_tcrs, full_mask))
            neighborhood_cdr3s = [s.split(',')[1] for s in neighborhood_tcrs]
            annulus_enrichments = list(compress(self.scores, annulus_mask))
            self.all_radii_dict[radius] = {
                "mean_enrichment": np.mean(neighborhood_enrichments),
                "annulus_enrichment": np.mean(annulus_enrichments) if annulus_enrichments != [] else None,
                "tcrs": neighborhood_tcrs,
                "cluster_size": int(full_mask.sum())
            }
            print("{}, {}".format(neighborhood_mask.sum(), np.mean(neighborhood_enrichments)))
        
            records = [SeqRecord(Seq(cdr3), id=tcr) for tcr, cdr3 in zip(neighborhood_tcrs, neighborhood_cdr3s)]

            alignment_infile = os.path.join(DIRECTORIES[TMP_OUTPUT], "cluster_cdr3s.fasta")

            with open(alignment_infile, "w") as output_handle:
                SeqIO.write(records, output_handle, "fasta")
            
            #hmmer_manager.run_hmmalign(alignment_infile)
            #hmmer_manager.run_hmmbuild()
            #hmmer_manager.run_hmmstat()

            #self.all_radii_dict[radius]['entropy'] = hmmer_manager.hmm_stats['relent']

            previous_radius = radius
        
        self.df = pd.DataFrame.from_dict(
            {
                radius: {
                    'cluster_size': v['cluster_size'],
                    'mean_enrichment': v['mean_enrichment'],
                    'annulus_enrichment': v['annulus_enrichment'],
                } for radius, v in self.all_radii_dict.items()
            }
        )
        self.df = self.df.transpose()
        self.df['radius'] = self.all_radii_dict.keys()

        self.df.loc[:, ('radius', 'annulus_enrichment')].to_csv(self.seg_csv_file, index=False)

        os.system(f'Rscript R/segmented_regression.R {outdir}')
        with open(os.path.join(outdir, 'breakpoint.txt')) as f:
            bp = [float(line.strip()) for line in f]
            if len(bp) != 1:
                raise Exception("Exactly one breakpoint was not found")
            else:
                self.breakpoint_radius = bp[0]


        self.cluster_dict = self.all_radii_dict[self.breakpoint_radius]
        self.is_in_cluster = [tcr in self.cluster_dict['tcrs'] for tcr in self.unique_tcrs]
        self.cluster_df = pd.DataFrame({'tcr': self.unique_tcrs, 'score': self.scores, 'is_in_cluster': self.is_in_cluster})
