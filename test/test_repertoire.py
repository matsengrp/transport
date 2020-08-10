import unittest

from python.repertoire import Repertoire

class TestRepertoire(unittest.TestCase):
    def test_repertoire(self):
        rep_1 = Repertoire(
            filename="test/data/test_seqs.tcrs",
            distribution_type="uniform",
            compute_distance_matrix=True,
            species="human"
        )
    
        self.assertEqual(rep_1.total_N, 5)
        self.assertEqual(rep_1.unique_N, 4)
    
