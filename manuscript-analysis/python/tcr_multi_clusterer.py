import os
import sys

import numpy as np

sys.path.append(os.getcwd())
from common.params import DIRECTORIES, TMP_OUTPUT
from python.hmmer_manager import HMMerManager
from python.tcr_clusterer import TCRClusterer
from python.tcr_scorer import TCRScorer

class TCRMultiClusterer():
    def __init__(self, file_1, file_2, species, outdir, max_cluster_count=10):
        self.species = species

        all_clusters = []
        
        initial_scorer = TCRScorer(file_1=file_1, file_2=file_2, species=self.species)
        initial_rep_2_self_dist_mat = initial_scorer.repertoire_2.distance_matrix
        initial_score_dict = initial_scorer.enrichment_dict
        cluster_outdir = os.path.join(outdir, "cluster_1")
        if not os.path.exists(cluster_outdir):
            os.makedirs(cluster_outdir)
        initial_clusterer = TCRClusterer(self_distance_matrix=initial_rep_2_self_dist_mat, score_dict=initial_score_dict, cluster_label="cluster_1", outdir=cluster_outdir)
        result = {tcr: {'score': score, 'cluster': 0} for tcr, score in initial_scorer.enrichment_dict.items()} 
        for tcr in initial_clusterer.cluster_dict['tcrs']:
            result[tcr]['cluster'] = 1
        
        sub_repertoire_tcrs = [tcr for tcr in initial_scorer.repertoire_2.unique_tcrs if tcr not in initial_clusterer.cluster_dict['tcrs']]
        hmmer_manager = HMMerManager()
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        hmmer_manager.build_hmm_from_sequences(
            [s.split(',')[1] for s in initial_clusterer.cluster_dict['tcrs']],
            outdir=cluster_outdir
        )
        cluster = 2
        sub_file = os.path.join(DIRECTORIES[TMP_OUTPUT], 'sub_rep.csv')
        while cluster <= max_cluster_count:
            np.savetxt(sub_file, sub_repertoire_tcrs, fmt="%s")
            score_dict = {k: initial_score_dict[k] for k in sub_repertoire_tcrs}
            cluster_outdir = os.path.join(outdir, f"cluster_{str(cluster)}")
            if not os.path.exists(cluster_outdir):
                os.makedirs(cluster_outdir)
            current_scorer, current_clusterer = self.run_clustering_step(file_1=file_1, file_2=sub_file, cluster_label=cluster, score_dict=score_dict, cluster_outdir=cluster_outdir) 
            sub_repertoire_tcrs = [tcr for tcr in current_scorer.repertoire_2.unique_tcrs if tcr not in current_clusterer.cluster_dict['tcrs']]
            for tcr in current_clusterer.cluster_dict['tcrs']:
                result[tcr]['cluster'] = cluster
            hmmer_manager = HMMerManager()
            print(cluster_outdir)
            hmmer_manager.build_hmm_from_sequences(
                [s.split(',')[1] for s in current_clusterer.cluster_dict['tcrs']],
                outdir=cluster_outdir
            )
            cluster = cluster + 1

        self.result = result

    def run_clustering_step(self, file_1, file_2, cluster_label, score_dict=None, cluster_outdir=None):
        scorer = TCRScorer(file_1=file_1, file_2=file_2, species=self.species)
        if score_dict is None:
            score_dict=scorer.enrichment_dict
        rep_2_self_dist_mat = scorer.repertoire_2.distance_matrix
        tcr_clusterer = TCRClusterer(self_distance_matrix=rep_2_self_dist_mat, score_dict=score_dict, cluster_label=f"cluster_{cluster_label}", outdir=cluster_outdir)
        return scorer, tcr_clusterer

