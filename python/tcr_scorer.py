import numpy as np
import ot

from common.params import DEFAULT_LAMBDA, DMAX
from python.repertoire import Repertoire
from python.tcr_dist import TCRDist

class TCRScorer():
    def __init__(self, file_1, file_2, distribution_type="uniform", lambd=DEFAULT_LAMBDA):
        self.file_1 = file_1
        self.file_2 = file_2
        self.distribution_type = distribution_type
        self.lambd = lambd

        self.repertoire_1 = Repertoire(self.file_1, self.distribution_type)
        self.repertoire_2 = Repertoire(self.file_2, self.distribution_type)
    
        self.dist_mat = TCRDist().get_raw_distance_matrix(
            self.repertoire_1.deduplicated_filename,
            self.repertoire_2.deduplicated_filename,
            verbose=False
        )/DMAX

        self.compute_efforts()
    
    def compute_efforts(self):
        self.ot_mat = ot.sinkhorn(self.repertoire_1.mass, self.repertoire_2.mass, self.dist_mat, self.lambd)
        self.effort_mat = np.multiply(self.dist_mat, self.ot_mat)
    
        N2 = self.repertoire_2.total_N
        efforts = DMAX*N2*self.effort_mat.sum(axis=0)
    
        assert len(efforts) == self.repertoire_2.unique_N
    
        self.effort_dict = {tcr: effort for tcr, effort in zip(self.repertoire_2.unique_tcrs, efforts)}
