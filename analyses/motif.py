import json
import os
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import pandas as pd

sys.path.append(os.getcwd())

from config import CONFIG
from python.hmmer_manager import HMMerManager
from python.tcr_clusterer import TCRClusterer
from python.utils import extract_cdr3s

with open(os.path.join(CONFIG["JSON_OUTPUT"], 'empirical_fg_bg_nbhd_stats.json')) as f:
    result = json.load(f)

subjects = result.keys()

sample_sizes = {
    subject: len(result[subject][str(CONFIG["DEFAULT_NEIGHBOR_RADIUS"])]) 
    for subject in subjects
}

sample_size_threshold = 1

def get_cluster_objects_from_subject(subject):
    subject_distance_matrix = np.loadtxt(os.path.join(CONFIG["DIST_MATRICES_OUTPUT"], subject + '.csv'), dtype='i', delimiter=',')
    info_dict = result[subject][str(CONFIG["DEFAULT_NEIGHBOR_RADIUS"])]
    score_dict = {tcr: tcr_info['foreground']['score'] for tcr, tcr_info in info_dict.items()}

    tcr_clusterer = TCRClusterer(subject_distance_matrix, score_dict)
    
    df = tcr_clusterer.df
    df['subject'] = subject

    return df, tcr_clusterer.all_radii_dict

def run_cluster_analysis():
    dfs = []
    full_dict = {}

    for subject in subjects:
        print(subject, sample_sizes[subject])
        if sample_sizes[subject] > sample_size_threshold:
            subject_cluster_df, subject_all_radii_dict = get_cluster_objects_from_subject(subject)
            dfs.append(subject_cluster_df)
            full_dict[subject] = subject_all_radii_dict

    full_df = pd.concat(dfs)
    full_df.to_csv(os.path.join(CONFIG["CSV_OUTPUT"], "motif.csv"), index=False)

    with open(os.path.join(CONFIG["JSON_OUTPUT"], "motif.json"), "w") as fp:
        json.dump(full_dict, fp)

def get_profile_from_subject_cluster(subject, cluster_radius):
    cluster_objects = get_cluster_objects_from_subject(subject + ".tcrs")
    cluster = cluster_objects[1][cluster_radius]['tcrs']
    cluster_cdr3s = extract_cdr3s(cluster) 

    cluster_cdr3s_fasta = os.path.join(CONFIG["TMP_OUTPUT"], subject + "_cluster_cdr3s.fasta")
    cluster_cdr3s_sto = os.path.join(CONFIG["TMP_OUTPUT"], subject + "_cluster_cdr3s.sto")
    cluster_cdr3s_hmm = os.path.join(CONFIG["TMP_OUTPUT"], subject + "_cluster_cdr3s.hmm")

    cluster_records = [SeqRecord(Seq(cdr3), id=cdr3) for cdr3 in cluster_cdr3s]
    with open(cluster_cdr3s_fasta, "w") as output_handle:
        SeqIO.write(cluster_records, output_handle, "fasta")

    hmmer_manager = HMMerManager()
    hmmer_manager.run_hmmalign(alignment_infile=cluster_cdr3s_fasta, alignment_outfile=cluster_cdr3s_sto, hmm_outfile=cluster_cdr3s_hmm)
    hmmer_manager.run_hmmbuild(hmm_file=cluster_cdr3s_hmm, alignment_outfile=cluster_cdr3s_sto)
    return cluster_cdr3s_hmm

def run_motif_analysis(reference_subject, cluster_radius, outfile_prefix):
    reference_cluster_cdr3s_hmm = get_profile_from_subject_cluster(reference_subject, cluster_radius) 

    for subject in subjects:
        if sample_sizes[subject] > sample_size_threshold:

            subject_all_cdr3s_fasta = os.path.join(CONFIG["TMP_OUTPUT"], "{}_all_cdr3s.fasta".format(subject))
            subject_hmmsearch_outfile = os.path.join(
                CONFIG["TMP_OUTPUT"], 
                "{}_hmmsearch.out".format(subject)
            )
            subject_df = pd.read_csv(
                os.path.join(CONFIG["IEL_DATA_DIR"], subject), 
                header=None, names=['v_gene', 'cdr3']
            )
            subject_all_cdr3s = subject_df['cdr3']
            subject_all_cdr3s_fasta = os.path.join(CONFIG["TMP_OUTPUT"], "{}_all_cdr3s.fasta".format(subject))
            subject_all_records = [SeqRecord(Seq(cdr3), id=cdr3) for cdr3 in subject_all_cdr3s]
            with open(subject_all_cdr3s_fasta, "w") as output_handle:
                SeqIO.write(subject_all_records, output_handle, "fasta")

            subject_hmmsearch_result = HMMerManager().run_hmmsearch(reference_cluster_cdr3s_hmm, subject_all_cdr3s, subject_hmmsearch_outfile, sequence_ids=subject_all_cdr3s)
            subject_e_values = [{'cdr3': tcr_info['target_name'], 'e_value': tcr_info['e_value']} for tcr_info in subject_hmmsearch_result]
            subject_e_value_df = pd.DataFrame(subject_e_values)
            subject_e_value_df['subject'] = subject
            e_value_dfs.append(subject_e_value_df)
    
    e_value_df = pd.concat(e_value_dfs)
    
    e_value_df.to_csv(os.path.join(CONFIG["CSV_OUTPUT"], "{}_e_value.csv".format(outfile_prefix)), index=False)

if __name__ == "__main__":
    if not os.path.exists(CONFIG["CSV_OUTPUT"]):
        os.makedirs(CONFIG["CSV_OUTPUT"])

    run_cluster_analysis()
