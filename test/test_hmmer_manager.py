import os
import sys
import unittest

from Bio import SeqIO
import numpy as np

sys.path.append(os.getcwd())

from python.hmmer_manager import HMMerManager

class TestHMMerManager(unittest.TestCase):
    def test_hmmer(self):
        hmmer_manager = HMMerManager()
        sequences = []
        with open("test/data/test_cdr3s.fasta", "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences.append(record.seq._data)

        outdir = "test/tmp_output"
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        hmmer_manager.build_hmm_from_sequences(sequences, outdir=outdir)
        result = hmmer_manager.run_hmmsearch(
            os.path.join(outdir, 'seqs.hmm'),
            os.path.join(outdir, 'aligned_seqs.sto'),
            os.path.join(outdir, 'test.out')
        )
        self.assertEqual(len(result), 6)
        self.assertTrue(np.all([float(i['e_value']) < 1e-08 for i in result]))

        # Check for an exception if the query sequence file does not exist
        with self.assertRaises(Exception):
            result = hmmer_manager.run_hmmsearch(
                os.path.join(outdir, 'seqs.hmm'),
                os.path.join(outdir, 'blah.sto'),
                os.path.join(outdir, 'test.out')
            )

        # Check for an exception if the hmmer file does not exist
        with self.assertRaises(Exception):
            result = hmmer_manager.run_hmmsearch(
                os.path.join(outdir, 'blah.hmm'),
                os.path.join(outdir, 'aligned_seqs.sto'),
                os.path.join(outdir, 'test.out')
            )
