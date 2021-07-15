import os

import numpy as np
import pandas as pd

class TCRDist():

    def __init__(self, species, species_db, tcrdists_exe):
        self.db = species_db
        self.exe = tcrdists_exe
        self.species = species
        #if self.species == "mouse":
        #    self.db = "pot_data/fake_pubtcrs_db_mouse"
        #elif self.species == "human":
        #    self.db = "pot_data/db"
        #else:
        #    raise Exception("Unsupported species (can only be mouse or human)")

    def get_raw_distance_matrix( 
        self,
        f1,
        f2,
        as_pandas_dataframe=False,
        index_column=None,
        verbose=True,
        output_dir="tmp_output",
    ):
    
        cmd = '{} -i {} -j {} -d {} --terse'.format(
        # cmd = '{} -i {} -j {} -d {} -g {}'.format(
            self.exe,
            os.path.join(output_dir, f1),
            os.path.join(output_dir, f2),
            self.db,
            # self.species
        )
        if verbose:
            print(cmd)
        all_dists = []
        for line in os.popen(cmd):
            
            try:
                all_dists.append( [float(x) for x in line.split() ] )
            except ValueError:
                print(line)
        #print(all_dists)
        #sys.exit()

        N1 = len(all_dists)
        N2 = len(all_dists[0])
        for dists in all_dists:
            assert len(dists) == N2
    
        D = np.array(all_dists)
        if verbose:
            print('loaded dists',D.shape)
    
        if as_pandas_dataframe:
            if index_column:
                D = pd.DataFrame(D)
    
                df_1 = get_df_from_file(f1)
                df_2 = get_df_from_file(f2)
    
                D.index = df_1.iloc[:, index_column]
                D.columns = df_2.iloc[:, index_column]
    
        return D
