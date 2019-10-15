import numpy as np
import os
import pandas as pd
import random

def compute_runtimes(dist_mat_filename):
    pass

if __name__ == "__main__":
    dist_mat_filename = "Data/all_IELrep_beta_distances.txt"
    batch_filename = "Data/tmp_distances.tsv"
    batch_sizes = [int(np.exp(x)) for x in range(3, 9)]
    for batch_size in batch_sizes:
        # First, subsample to batch_size rows using shuf, avoiding need to load full matrix into memory
        command = "shuf -n " + str(batch_size) + " " + dist_mat_filename + " > " + batch_filename
        os.system(command)
        dist_mat = pd.read_csv(batch_filename, sep=" ", header=None)
        num_columns = dist_mat.shape[1]
        # Next, subsample to batch_size columns
        dist_mat = dist_mat.iloc[:, random.sample(range(0, num_columns), batch_size)]
