import numpy as np
import os
import ot
import pandas as pd
import random
import time

def compute_runtimes(dist_mat_filename):
    pass

if __name__ == "__main__":
    dist_mat_filename = "Data/all_IELrep_beta_distances.txt"
    batch_filename = "Data/tmp_distances.tsv"
    batch_sizes = [int(np.exp(x)) for x in range(3, 11)]
    for batch_size in batch_sizes:
        # First, subsample to batch_size rows using shuf, avoiding need to load full matrix into memory
        command = "shuf -n " + str(batch_size) + " " + dist_mat_filename + " > " + batch_filename
        os.system(command)
        dist_mat = pd.read_csv(batch_filename, sep=" ", header=None)
        num_columns = dist_mat.shape[1]
        # Next, subsample to batch_size columns
        dist_mat = dist_mat.iloc[:, random.sample(range(0, num_columns), batch_size)]

        start_time = time.time()
        lambdas = np.exp(np.linspace(-1 ,3, 20))
        empirical_distribution = np.ones((batch_size,)) / batch_size
        dists = [ot.sinkhorn2(empirical_distribution, empirical_distribution, dist_mat, lam) for lam in lambdas]
        mats = [ot.sinkhorn(empirical_distribution, empirical_distribution, dist_mat, lam) for lam in lambdas]
        elapsed_time = time.time() - start_time
        print(str(batch_size) + ": " + str(elapsed_time))
