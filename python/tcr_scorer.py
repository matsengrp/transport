from collections import defaultdict

import numpy as np
import ot

from common.params import DEFAULT_LAMBDA, DEFAULT_NEIGHBOR_RADIUS, DMAX
from python.repertoire import Repertoire
from python.tcr_dist import TCRDist

class TCRScorer():
    def __init__(self, file_1, file_2, distribution_type="uniform", lambd=DEFAULT_LAMBDA, neighbor_radius=DEFAULT_NEIGHBOR_RADIUS):
        self.file_1 = file_1
        self.file_2 = file_2
        self.distribution_type = distribution_type
        self.lambd = lambd
        self.neighbor_radius = neighbor_radius

        self.repertoire_1 = Repertoire(self.file_1, self.distribution_type)
        self.repertoire_2 = Repertoire(self.file_2, self.distribution_type, compute_distance_matrix=True)
    
        self.dist_mat = TCRDist().get_raw_distance_matrix(
            self.repertoire_1.deduplicated_filename,
            self.repertoire_2.deduplicated_filename,
            verbose=False
        )/DMAX

        self.compute_efforts()
        self.compute_enrichments()
    
    def compute_efforts(self):
        self.ot_mat = ot.sinkhorn(self.repertoire_1.mass, self.repertoire_2.mass, self.dist_mat, self.lambd)
        self.effort_matrix = np.multiply(self.dist_mat, self.ot_mat)
        N2 = self.repertoire_2.total_N
        self.efforts = DMAX*N2*self.effort_matrix.sum(axis=0)
        assert len(self.efforts) == self.repertoire_2.unique_N
    
        self.effort_dict = {tcr: effort for tcr, effort in zip(self.repertoire_2.unique_tcrs, self.efforts)}

    def compute_enrichments(self, neighbor_radius=DEFAULT_NEIGHBOR_RADIUS):
        self.enrichment_dict = defaultdict()
        self.neighbor_counts = defaultdict()
        for i, tcr in zip(range(self.repertoire_2.unique_N), self.repertoire_2.unique_tcrs):
            neighborhood_mask = ( self.repertoire_2.distance_matrix[i, :] < neighbor_radius )
            self.enrichment_dict[tcr] = self.efforts[neighborhood_mask].sum()
            self.neighbor_counts[tcr] = neighborhood_mask.sum()
