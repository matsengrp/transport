from glob import glob
import os
import sys

import pandas as pd

sys.path.append(os.getcwd())
from common.params import DIRECTORIES, IEL_DATA_DIR, TMP_OUTPUT
from python.tcr_scorer import TCRScorer
from python.utils import get_df_from_file

def write_full_replicate_dataset(output_filename, reptag, chain="B"):
    repfiles = sorted(glob('{}{}_*_{}.tcrs'.format(IEL_DATA_DIR, reptag, chain)))
    dfs = []
    for repfile in repfiles:
        dfs.append(get_df_from_file(repfile))
    full_df = pd.concat(dfs)
    full_df.iloc[:, [0, 1]].to_csv(output_filename, index=False, header=False)

cd4_reptag = "CD4"
dn_reptag = "DN"
chain = "B"

dn_file = os.path.join(DIRECTORIES[TMP_OUTPUT], "all_dn.csv")
cd4_file = os.path.join(DIRECTORIES[TMP_OUTPUT], "all_cd4.csv")
write_full_replicate_dataset(dn_file, "DN")
write_full_replicate_dataset(cd4_file, "CD4")

scorer = TCRScorer(file_1=cd4_file, file_2=dn_file)
import pdb; pdb.set_trace()
