from glob import glob
import json
import os
import subprocess
import sys

import pandas as pd

sys.path.append(os.getcwd())
from common.params import CSV_OUTPUT_DIRNAME, DIRECTORIES, HMM_OUTPUT, IEL_DATA_DIR, JSON_OUTPUT, TMP_OUTPUT
from python.hmmer_manager import HMMerManager
from python.tcr_clusterer import TCRClusterer
from python.tcr_dist import TCRDist
from python.tcr_multi_clusterer import TCRMultiClusterer
from python.tcr_scorer import TCRScorer
from python.utils import get_df_from_file

bulk_a = "../../pot_data/spleen-data/bulk_mouse_spleen_la_mc_pcr_A.csv"
bulk_b = "../../pot_data/spleen-data/bulk_mouse_spleen_la_mc_pcr_B.csv"

multi_clusterer = TCRMultiClusterer(file_1=bulk_a, file_2=bulk_b, species="mouse", outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "spleen"))

pd.DataFrame(multi_clusterer.result).transpose().to_csv(os.path.join(DIRECTORIES[TMP_OUTPUT], 'multiple_clusters_spleen.csv'))
