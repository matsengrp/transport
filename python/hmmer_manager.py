from io import BytesIO
import json
import os

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from PIL import Image
import requests

from common.params import DIRECTORIES, TMP_OUTPUT, TRB_MOUSE_CDR3_HMM, TRB_MOUSE_CDR3_STO

class HMMerManager():
    hmm_filename = TRB_MOUSE_CDR3_HMM
    base_alignment_filename = TRB_MOUSE_CDR3_STO

    if not os.path.exists(hmm_filename):
        os.mknod(hmm_filename)

    def __init__(self):
        self.alignment_outfile = os.path.join(DIRECTORIES[TMP_OUTPUT], "cluster_cdr3s.sto")
        self.motif_hmm_file = os.path.join(DIRECTORIES[TMP_OUTPUT], "motif.hmm")
        self.motif_hmm_stats_file = os.path.join(DIRECTORIES[TMP_OUTPUT], "motif_stats.txt")

    def run_hmmalign(self, alignment_infile, alignment_outfile=None):
        if alignment_outfile is None:
            alignment_outfile = self.alignment_outfile

        command = 'hmmalign {} {} > {}'.format(
            self.hmm_filename,
            alignment_infile,
            alignment_outfile
        )
        print(command)
        os.system(command)

    def run_hmmbuild(self, hmm_file=None, alignment_outfile=None):
        if hmm_file is None:
            hmm_file = self.motif_hmm_file
        else:
            self.motif_hmm_file = hmm_file

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


    def build_hmm_from_sequences(
            self,
            sequence_list,
            hmm_filename=None,
            alignment_outfilename=None,
            fasta_filename=None,
            plot_filename=None
        ):
        if fasta_filename is None:
            fasta_filename = os.path.join(DIRECTORIES[TMP_OUTPUT], 'seqs.fasta')

        if alignment_outfilename is None:
            alignment_outfilename = os.path.join(DIRECTORIES[TMP_OUTPUT], 'seqs.sto')

        if hmm_filename is None:
            hmm_filename = os.path.join(DIRECTORIES[TMP_OUTPUT], 'seqs.hmm')

        if plot_filename is None:
            plot_filename = os.path.join(DIRECTORIES[TMP_OUTPUT], 'logo_plot.png')
                
        records = [SeqRecord(Seq(sequence), id=sequence) for sequence in sequence_list]
        with open(fasta_filename, 'w') as output_handle:
            SeqIO.write(records, output_handle, "fasta")

        self.run_hmmalign(alignment_infile=fasta_filename, alignment_outfile=alignment_outfilename)
        self.run_hmmbuild(hmm_file=hmm_filename, alignment_outfile=alignment_outfilename)
        self.get_logo_plot(plot_outfilename=plot_filename)

    def get_logo_plot(self, plot_outfilename=None):
        if plot_outfilename is None:
            plot_outfilename = os.path.join(DIRECTORIES[TMP_OUTPUT], 'logo_plot.png')

        r = requests.post(
            'http://skylign.org',
            headers={"Accept": "application/json"},
            files={"file": open(self.motif_hmm_file, "rb")},
            data={"processing": "hmm"}
        )
        r2 = requests.get(
            json.loads(r.text)['url'],
            headers={"Accept": "image/png"}
        )
        i = Image.open(BytesIO(r2.content))
        i.save(plot_outfilename)

