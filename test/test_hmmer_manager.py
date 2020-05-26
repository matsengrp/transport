import os
import sys

import numpy as np

sys.path.append(os.getcwd())

from python.hmmer_manager import HMMerManager

def test_hmmer():
    hmmer_manager = HMMerManager()
    hmmer_manager.run_hmmalign('test/data/test_cdr3s.fasta', alignment_outfile='tmp_output/test.sto')
    hmmer_manager.run_hmmbuild('tmp_output/test.hmm', 'tmp_output/test.sto')
    result = hmmer_manager.run_hmmsearch('tmp_output/test.hmm', 'tmp_output/test.sto', 'tmp_output/test.out')
    assert len(result) == 6
    assert np.all([float(i['e_value']) < 1e-08 for i in result])
