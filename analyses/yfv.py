from collections import defaultdict
from glob import glob
import json
import os
import sys

import numpy as np
import ntpath

sys.path.append(".")
from common.params import DIRECTORIES, HMM_OUTPUT
from python.hmmer_manager import HMMerManager
from python.tcr_multi_clusterer import TCRMultiClusterer
from config import CONFIG

seq_data_dir = CONFIG["YFV_DATA_DIR_TOP1000"]

def run_yfv_analysis(file_1, file_2, outdir):
    tcr_multi_cluster = TCRMultiClusterer(file_1=file_1, file_2=file_2, species="human", outdir=outdir)

run_yfv_analysis(
    os.path.join(seq_data_dir, "P1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "P1_pre0_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "P1_0_-7")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "P1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "P1_0_F2_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "P1_0_0")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "P1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "P1_15_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "P1_0_15")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "P1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "P1_45_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "P1_0_45")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "P2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "P2_pre0_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "P2_0_-7")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "P2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "P2_0_F2_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "P2_0_0")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "P2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "P2_15_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "P2_0_15")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "P2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "P2_45_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "P2_0_45")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "Q1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "Q1_pre0_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "Q1_0_-7")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "Q1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "Q1_0_F2_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "Q1_0_0")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "Q1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "Q1_15_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "Q1_0_15")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "Q1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "Q1_45_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "Q1_0_45")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "Q2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "Q2_pre0_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "Q2_0_-7")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "Q2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "Q2_0_F2_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "Q2_0_0")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "Q2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "Q2_15_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "Q2_0_15")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "Q2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "Q2_45_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "Q2_0_45")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "S1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "S1_pre0_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "S1_0_-7")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "S1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "S1_0_F2_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "S1_0_0")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "S1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "S1_15_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "S1_0_15")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "S1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "S1_45_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "S1_0_45")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "S2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "S2_pre0_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "S2_0_-7")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "S2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "S2_0_F2_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "S2_0_0")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "S2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "S2_15_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "S2_0_15")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "S2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "S2_45_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(CONFIG["HMM_OUTPUT"], "S2_0_45")
)
