## need to fixup the variable names and output columns now that I've looked back at
##  the terminology in the paper

import numpy as np
import ot
from glob import glob
from os import popen, remove
from os.path import exists
import sys
import random
import pandas as pd

## this is the tcrdist-computing executable
exe = 'pubtcrs/bin/tcrdists'

## this is a db-directory needed for the tcrdists calc
db = 'data/human_db'

naive_tcrs_file = 'data/flu_spike_ins/naive_tcrs.tsv'
flu_tcrs_file = 'data/flu_spike_ins/flu_tcrs.tsv'

results_file = 'results/flu_spike_ins/flu_auc_values.tsv'

#lambd = .025
lambd = .01
random_state = 10
#nbrcutoff_single = 48.5
nbrcutoff_paired = 120.5
#Dmax_single = 200 # constant across comparisons
Dmax_paired = 400 # constant across comparisons

num_naive_list = [1000]
num_flu_list = [5, 10, 20, 40, 80, 160, 320]
num_repeats = 10

tmpfile_prefix = 'tmpfile_flu_spike_ins'

######################################################################################88
######################################################################################88
######################################################################################88

naive_tcrs_df = pd.read_table(naive_tcrs_file)
flu_tcrs_df = pd.read_table(flu_tcrs_file)


cleanup_tmpfiles = True

tmpfiles = []
def get_tmpfile():
    global tmpfiles
    tmpfile = f'{tmpfile_prefix}_{random.random()}.txt'
    tmpfiles.append(tmpfile) # will delete later, maybe
    return tmpfile

def compute_auroc(scores, labels, ascending=True):
    ''' lower scores are better (if ascending==True)
    labels are True for positive, False for negative

    heuristic tie handling with randomness

    '''
    sortl = sorted(zip(scores, np.random.rand(len(scores)), labels),
                   reverse=not ascending)

    auc = 0
    num_correct, num_incorrect = 0, 0
    total_correct = sum(x[2] for x in sortl)
    total_incorrect = sum(not x[2] for x in sortl)

    for score, _, label in sortl:
        if label:
            num_correct += 1
        else:
            num_incorrect += 1
            auc += num_correct/(total_correct*total_incorrect)
    return auc


def get_raw_distance_matrix( f1, f2 ):
    cmd = f'{exe} -i {f1} -j {f2} -d {db} --terse'
    #print(cmd)
    all_dists = []
    for line in popen(cmd):
        all_dists.append( [float(x) for x in line.split() ] )
    N1 = len(all_dists)
    N2 = len(all_dists[0])
    for dists in all_dists:
        assert len(dists) == N2

    D = np.array(all_dists)
    #print('loaded dists',D.shape)
    return D


def get_paired_distance_matrix(df1, df2):
    f1 = get_tmpfile()
    f2 = get_tmpfile()
    matrices = []
    for ab in 'ab':
        cols = f'v{ab} cdr3{ab}'.split()
        df1[cols].to_csv(f1, header=False, index=False)
        df2[cols].to_csv(f2, header=False, index=False)
        matrices.append(get_raw_distance_matrix(f1,f2))
    return matrices[0] + matrices[1]


dfl = []
for num_naive in num_naive_list:
    for num_flu in num_flu_list:
        for r in range(num_repeats):
            print(num_naive, num_flu)
            naive_sample = naive_tcrs_df.sample(
                2*num_naive-num_flu, random_state=random_state+r)
            fg_df = pd.concat([
                flu_tcrs_df.sample(num_flu, random_state=random_state+r),
                naive_sample.iloc[num_naive:],
            ])
            bg_df = naive_sample.iloc[:num_naive]
            D_fg_bg = get_paired_distance_matrix(fg_df, bg_df)
            D_fg_fg = get_paired_distance_matrix(fg_df, fg_df)

            N1 = num_naive # foreground
            N2 = num_naive # background (in case we change sizes later)

            assert D_fg_bg.shape == (N1,N2)
            assert D_fg_fg.shape == (N1,N1)

            D = np.array( D_fg_bg ) # make copy, then normalize
            D /= Dmax_paired
            wts1 = np.ones((N1,)) / N1
            wts2 = np.ones((N2,)) / N2

            print('run sinkhorn2 for dist',N1,N2)
            dist = ot.sinkhorn2(wts1, wts2, D, lambd )
            dist *= Dmax_paired

            # is this silly to run twice?? yup, dist is just (D*mat).sum()...
            print('run sinkhorn for mat',N1,N2)
            mat = ot.sinkhorn(wts1, wts2, D, lambd)

            # total effort matrix
            hadamard = Dmax_paired * np.multiply(D, mat)


            print(f'sinkhorn_dist: {dist:9.3f} {hadamard.sum():9.3f} '
                  f'lambda: {lambd} Dmax: {Dmax_paired} '
                  f' Ns: {N1} {N2}', flush=True)

            ## take row and column (unused) sums of the effort matrix
            ## this is TotalLoneliness (finally looking back at the paper now)
            fg_transport = hadamard.sum(axis=1)
            assert fg_transport.shape[0] == N1


            # look through the file1 tcrs
            # get total transport in nbrhood
            nbrs = D_fg_fg <= nbrcutoff_paired
            nbrhood_sizes = nbrs.sum(axis=-1)

            # this is RelativeLoneliness or just loneliness
            fg_loneliness = (nbrs * fg_transport[None,:]).sum(axis=1)

            labels = np.arange(N1) < num_flu

            # auc wrt to RelativeLoneliness aka loneliness
            auc_loneliness = compute_auroc(fg_loneliness, labels, ascending=False)

            # auc wrt to TotalLoneliness
            auc_transport = compute_auroc(fg_transport, labels, ascending=False)

            print(f'auc_loneliness: {auc_loneliness:.3f} '
                  f'auc_transport: {auc_transport:.3f} '
                  f'num_naive: {num_naive:5d} num_flu: {num_flu:5d}',
                  flush=True)

            dfl.append(dict(
                num_naive=num_naive,
                num_flu=num_flu,
                auc_loneliness=auc_loneliness,
                auc_transport=auc_transport,
                repeat=r,
            ))

            pd.DataFrame(dfl).to_csv(results_file, sep='\t', index=False)
            continue


            for ii in np.argsort(fg_loneliness)[:-50:-1]:
                row = fg_df.iloc[ii]
                print(f'f1_index: {ii:4d} transport: {fg_transport[ii]:9.3f} '
                      f'loneliness: {fg_loneliness[ii]:9.3f} '
                      f'num_nbrs:: {nbrhood_sizes[ii]-1:2d} '
                      f'tcr: {row.va} {row.cdr3a} {row.vb} {row.cdr3b}')

            for ii in np.argsort(fg_transport)[:-50:-1]:
                row = fg_df.iloc[ii]
                print(f'f1_index: {ii:4d} transport: {fg_transport[ii]:9.3f} '
                      f'loneliness: {fg_loneliness[ii]:9.3f} '
                      f'num_nbrs:: {nbrhood_sizes[ii]-1:2d} '
                      f'tcr: {row.va} {row.cdr3a} {row.vb} {row.cdr3b}')

pd.DataFrame(dfl).to_csv(results_file, sep='\t', index=False)

if cleanup_tmpfiles:
    for f in tmpfiles:
        if exists(f):
            remove(f)


