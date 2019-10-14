import numpy as np
import os
import pandas as pd

def compute_runtimes(dist_mat_filename):
    pass


if __name__ == "__main__":
    dist_mat_filename = "Data/all_IELrep_beta_distances.txt"
    batch_filename = "Data/tmp_distances.tsv"
    batch_sizes = [np.floor(np.exp(x)) for x in range(3, 9)]
    for batch_size in batch_sizes:
        command = "shuf -n " + str(batch_size) + " " + dist_mat_filename + " > " + batch_filename
        os.system(command)
        dist_mat = pd.read_csv(batch_filename, sep=" ", header=None)
