import os

from common.params import DIRECTORIES, TMP_OUTPUT, TRB_MOUSE_HMM

class HMMerManager():
    hmm_filename = TRB_MOUSE_HMM

    if not os.path.exists(hmm_filename):
        os.mknod(hmm_filename)

    def __init__(self):
        self.alignment_outfile = os.path.join(DIRECTORIES[TMP_OUTPUT], "cluster_cdr3s.sto")
        self.motif_hmm_file = os.path.join(DIRECTORIES[TMP_OUTPUT], "motif.hmm")
        self.motif_hmm_stats_file = os.path.join(DIRECTORIES[TMP_OUTPUT], "motif_stats.txt")

    def run_hmmalign(self, alignment_infile, alignment_outfile=None):
        if alignment_outfile is None:
            alignment_outfile = self.alignment_outfile

        command = 'hmmalign {} {} > {}'.format(self.hmm_filename, alignment_infile, alignment_outfile)
        print(command)
        os.system(command)

    def run_hmmbuild(self, hmm_file=None, alignment_outfile=None):
        if hmm_file is None:
            hmm_file = self.motif_hmm_file

        if alignment_outfile is None:
            alignment_outfile = self.alignment_outfile

        command = 'hmmbuild {} {}'.format(hmm_file, alignment_outfile)
        print(command)
        os.system(command)

    def run_hmmsearch(self, hmm_filename, sequence_database, outfile):
        command = 'hmmsearch --tblout {} -E 10000 {} {}'.format(outfile, hmm_filename, sequence_database)
        print(command)
        os.system(command)

        fields = ['target_name', 'accession', 'query_name', 'accession', 'e_value', 'score', 'bias', 'e_value_2', 'score_2', 'bias_2', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description_of_target']

        hmmsearch_result = []
        with open(outfile, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    values = line.split()
                    hmmsearch_result.append({field: value for field, value in zip(fields, values)})

        return hmmsearch_result

    def run_hmmstat(self):
        command = 'hmmstat {} > {}'.format(self.motif_hmm_file, self.motif_hmm_stats_file)
        os.system(command)

        fields = ['idx', 'name', 'accession', 'nseq', 'eff_nseq', 'M', 'relent', 'info', 'p relE', 'compKL']
        with open(self.motif_hmm_stats_file, 'r') as f:
            for line in f:
                print(line)
                if line.startswith('1'):
                    values = line.split()
                    self.hmm_stats = {field: value for field, value in zip(fields, values)}

