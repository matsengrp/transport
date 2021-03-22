from collections import Counter
import os
import re

from common.params import DIRECTORIES, TMP_OUTPUT, TCRDISTS_EXE, SPECIES_DB
from python.tcr_dist import TCRDist
from python.utils import get_df_from_file, write_deduplicated_file

class Repertoire():
    def __init__(self, filename, distribution_type, compute_distance_matrix=False, species=None):
        self.filename = filename
        self.distribution_type = distribution_type
        self.data_frame = get_df_from_file(self.filename)
        self.total_N = self.data_frame.shape[0]

        tcr_counter = Counter(self.data_frame['tcr'])
        self.unique_tcrs = tcr_counter.keys()
        self.unique_N = len(self.unique_tcrs)
        self.get_mass_distribution(tcr_counter)

        self.deduplicated_filename = re.sub(
            ".csv",
            "_deduplicated.csv", 
            re.sub(
                ".tcrs",
                "_deduplicated.tcrs",
                os.path.basename(self.filename)
            )
        )
        write_deduplicated_file(df=self.data_frame, filename=self.deduplicated_filename)

        if compute_distance_matrix:
            if species is None:
                raise Exception("species must be supplied if compute_distance_matrix is True")
            db = SPECIES_DB[species]
            self.distance_matrix = TCRDist(species=species, species_db=db, tcrdists_exe=TCRDISTS_EXE).get_raw_distance_matrix(self.deduplicated_filename, self.deduplicated_filename)

    def get_mass_distribution(self, tcr_counter):
        if self.distribution_type == "inverse_to_v_gene":
            def get_gene_masses(gene_list):
                unique_genes = list(set(gene_list))
                gene_mass_dict = {gene: 1/len(unique_genes) for gene in unique_genes}
                return gene_mass_dict
    
            self.mass = get_gene_weighted_mass_distribution(df)
            self.gene_mass_dict = get_gene_masses(df['v_gene'])
            return (mass, gene_mass_dict)
        elif self.distribution_type == "uniform":
            N = self.total_N
            self.mass = [count/N for count in tcr_counter.values()]
        else:
            raise Exception("Unsupported distribution_type")
