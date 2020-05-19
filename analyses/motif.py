import json
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pandas as pd

sys.path.append(os.getcwd())

from common.params import CSV_OUTPUT_DIRNAME, DEFAULT_NEIGHBOR_RADIUS, DIRECTORIES, DIST_MATRICES, JSON_OUTPUT, TMP_OUTPUT
from python.tcr_clusterer import HMMerManager, TCRClusterer

def get_cluster_objects_from_subject(subject):
    subject_distance_matrix = np.loadtxt(os.path.join(DIRECTORIES[DIST_MATRICES], subject + '.csv'), dtype='i', delimiter=',')
    info_dict = result[subject][str(DEFAULT_NEIGHBOR_RADIUS)]
    score_dict = {tcr: tcr_info['foreground']['score'] for tcr, tcr_info in info_dict.items()}

    tcr_clusterer = TCRClusterer(subject_distance_matrix, score_dict)
    
    df = tcr_clusterer.df
    df['subject'] = subject

    return df, tcr_clusterer.motif_dict

if __name__ == "__main__":
    with open(os.path.join(DIRECTORIES[JSON_OUTPUT], 'empirical_fg_bg_nbhd_stats.json')) as f:
        result = json.load(f)
    
    subjects = result.keys()
    sample_sizes = {subject: len(result[subject][str(DEFAULT_NEIGHBOR_RADIUS)]) for subject in subjects}
    sample_size_threshold = 300
    dfs = []
    full_dict = {}

    dn_12_result = get_cluster_objects_from_subject("DN_12_B.tcrs")
    dn_12_cluster = dn_12_result[1][70.5]['tcrs'] # Radius obtained from breakpoint script in R
    dn_12_cluster_cdr3s = [s.split(',')[1] for s in dn_12_cluster]
    dn_12_all_cdr3s = [s.split(',')[1] for s in list(result["DN_12_B.tcrs"][str(70.5)].keys())]

    dn_12_cluster_cdr3s_fasta = os.path.join(DIRECTORIES[TMP_OUTPUT], "dn_12_cluster_cdr3s.fasta")
    dn_12_cluster_cdr3s_sto = os.path.join(DIRECTORIES[TMP_OUTPUT], "dn_12_cluster_cdr3s.sto")
    dn_12_cluster_cdr3s_hmm = os.path.join(DIRECTORIES[TMP_OUTPUT], "dn_12_cluster_cdr3s.hmm")

    dn_12_all_cdr3s_fasta = os.path.join(DIRECTORIES[TMP_OUTPUT], "dn_12_all_cdr3s.fasta")
    dn_12_all_cdr3s_sto = os.path.join(DIRECTORIES[TMP_OUTPUT], "dn_12_all_cdr3s.sto")
    dn_12_all_cdr3s_hmm = os.path.join(DIRECTORIES[TMP_OUTPUT], "dn_12_all_cdr3s.hmm")

    dn_12_hmmsearch_outfile = os.path.join(DIRECTORIES[TMP_OUTPUT], "hmmsearch.out")

    all_records = [SeqRecord(Seq(cdr3), id=cdr3) for cdr3 in dn_12_all_cdr3s]
    with open(dn_12_all_cdr3s_fasta, "w") as output_handle:
        SeqIO.write(all_records, output_handle, "fasta")

    cluster_records = [SeqRecord(Seq(cdr3), id=cdr3) for cdr3 in dn_12_cluster_cdr3s]
    with open(dn_12_cluster_cdr3s_fasta, "w") as output_handle:
        SeqIO.write(cluster_records, output_handle, "fasta")

    dn_12_hmmer = HMMerManager()
    dn_12_hmmer.run_hmmalign(alignment_infile=dn_12_cluster_cdr3s_fasta, alignment_outfile=dn_12_cluster_cdr3s_sto)
    dn_12_hmmer.run_hmmbuild(hmm_file=dn_12_cluster_cdr3s_hmm, alignment_outfile=dn_12_cluster_cdr3s_sto)
    hmmsearch_result = dn_12_hmmer.run_hmmsearch(dn_12_cluster_cdr3s_hmm, dn_12_all_cdr3s_fasta, dn_12_hmmsearch_outfile)

    for subject in subjects:
        if sample_sizes[subject] > sample_size_threshold:
            subject_cluster_df, subject_motif_dict = get_cluster_objects_from_subject(subject)
            dfs.append(subject_cluster_df)
            full_dict[subject] = subject_motif_dict
    
    full_df = pd.concat(dfs)
    
    if not os.path.exists(CSV_OUTPUT_DIRNAME):
        os.makedirs(CSV_OUTPUT_DIRNAME)
    
    full_df.to_csv(os.path.join(CSV_OUTPUT_DIRNAME, "motif.csv"), index=False)
    
    with open(os.path.join(DIRECTORIES[JSON_OUTPUT], "motif.json"), "w") as fp:
        json.dump(full_dict, fp)
