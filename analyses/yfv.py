from collections import defaultdict
import json

import numpy as np
import ntpath

from glob import glob
import os
from scipy.stats import mannwhitneyu, ttest_ind
import sys

sys.path.append(os.getcwd())
from common.params import DIRECTORIES, HMM_OUTPUT
from python.hmmer_manager import HMMerManager
from python.tcr_multi_clusterer import TCRMultiClusterer
from python.tcr_clusterer import TCRClusterer
from python.tcr_scorer import TCRScorer

seq_data_dir = 'data/yfv'

def run_yfv_analysis(file_1, file_2, outdir):
    tcr_multi_cluster = TCRMultiClusterer(file_1=file_1, file_2=file_2, species="human", outdir=outdir)

run_yfv_analysis(
    os.path.join(seq_data_dir, "P1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "P1_0_F2_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "P1_0_0")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "P1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "P1_15_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "P1_0_15")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "P1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "P1_45_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "P1_0_45")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "P2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "P2_0_F2_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "P2_0_0")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "P2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "P2_15_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "P2_0_15")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "P2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "P2_45_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "P2_0_45")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "Q1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "Q1_0_F2_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "Q1_0_0")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "Q1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "Q1_15_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "Q1_0_15")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "Q1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "Q1_45_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "Q1_0_45")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "Q2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "Q2_0_F2_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "Q2_0_0")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "Q2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "Q2_15_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "Q2_0_15")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "Q2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "Q2_45_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "Q2_0_45")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "S1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "S1_0_F2_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "S1_0_0")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "S1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "S1_15_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "S1_0_15")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "S1_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "S1_45_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "S1_0_45")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "S2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "S2_0_F2_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "S2_0_0")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "S2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "S2_15_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "S2_0_15")
)
run_yfv_analysis(
    os.path.join(seq_data_dir, "S2_0_F1_.txt.top1000.tcrs"),
    os.path.join(seq_data_dir, "S2_45_F1_.txt.top1000.tcrs"),
    outdir=os.path.join(DIRECTORIES[HMM_OUTPUT], "S2_0_45")
)
