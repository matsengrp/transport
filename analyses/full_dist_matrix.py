import json

import numpy as np
import pandas as pd

from glob import glob
import os
import sys

sys.path.append(os.getcwd())
from common.params import DIRECTORIES, DIST_MATRICES
from python.utils import get_df_from_file, get_raw_distance_matrix, write_deduplicated_file

## this is the tcrdist-computing executable
exe = 'bin/tcrdists' 

## this is a db-directory needed for the tcrdists calc, for mouse tcrs
db = '/fh/fast/matsen_e/bolson2/transport/iel_data/fake_pubtcrs_db_mouse'


# on rhino1:
seq_data_dir = '/loc/no-backup/pbradley/share/pot_data/iels_tcrs_by_mouse/'

bg_reptag = 'DN'
chain = 'B'

bg_repfiles = sorted(glob('{}{}_*_{}.tcrs'.format(seq_data_dir, bg_reptag, chain)))

dfs = []
for repfile1 in bg_repfiles:
    dfs.append(get_df_from_file(repfile1))
    
full_df = pd.concat(dfs)

## compute intra-repertoire distance matrix to find TCR neighborhoods
dedup_repfile1 = "repfile_dedup.csv"
output_dir = "tmp_output"
write_deduplicated_file(full_df, dedup_repfile1, output_dir=output_dir)
D_11 = get_raw_distance_matrix(dedup_repfile1, dedup_repfile1)
np.savetxt(
    os.path.join(
        DIRECTORIES[DIST_MATRICES],
        "all_subjects.csv", 
    ),
    D_11,
    delimiter=",",
    fmt='%i'
)
