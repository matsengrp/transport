## this script compares per-tcr-efforts (ie row sums of the effort matrix) within a set of 'foreground'
## repertoires and between the foreground repertoires and a set of background repertoires. Here the
## foreground repertoires are the DN beta repertoires across the mice and the background reps are the
## CD4s (see fg_reptag and bg_reptag below)
##
## The idea is that each foreground tcr will have a distribution of efforts in the comparisons
## to the other foreground repertoires and also a distribution of efforts in the comparisons to
## the background repertoires. For a wonky tcr that's an outlier, both those distributions will
## be large. For a tcr that has more neighbors in the foreground repertoires than in the background
## repertoires, the two distributions will be different and we can attach a significance to that.
##
## Finally we are comparing the straight tcr-efforts with tcr-efforts summed over neighborhoods,
## to see which shows more significant differences.
##
##


from collections import defaultdict
import json

import numpy as np
import ntpath

import ot
from glob import glob
from os import popen
import sys
from scipy.stats import mannwhitneyu, ttest_ind


## this is the tcrdist-computing executable
exe = 'bin/tcrdists' #'/home/pbradley/gitrepos/pubtcrs/bin/tcrdists'
#exe = 'bin/tcrdists'

## this is a db-directory needed for the tcrdists calc, for mouse tcrs
db = '/fh/fast/matsen_e/bolson2/transport/iel_data/fake_pubtcrs_db_mouse' #'/loc/no-backup/pbradley/share/pot_data/fake_pubtcrs_db_mouse'

#db = 'data/databases/fake_pubtcrs_db_mouse'


# on rhino1:
seq_data_dir = '/loc/no-backup/pbradley/share/pot_data/iels_tcrs_by_mouse/'


## fg = foreground, bg = background
fg_reptag = 'DN'
bg_reptag = 'CD4'
chain = 'B'

fg_repfiles = sorted(glob('{}{}_*_{}.tcrs'.format(seq_data_dir, fg_reptag, chain)))
bg_repfiles = sorted(glob('{}{}_*_{}.tcrs'.format(seq_data_dir, bg_reptag, chain)))

lambd = .01
Dmax = 200 # constant across comparisons


nbrcutoffs = [i + .5 for i in range(0, 100, 5)] # distance at which two single-chain tcrs are considered nbrs



def get_raw_distance_matrix( f1, f2 ):
    ''' Returns the tcrdist distance matrix D between tcrs in file f1 and tcrs in file f2
    D.shape = (num_f1_tcrs, num_f2_tcrs)
    '''
    cmd = '{} -i {} -j {} -d {} --terse'.format( exe, f1, f2, db )
    print(cmd)
    all_dists = []
    for line in popen(cmd):
        all_dists.append( [float(x) for x in line.split() ] )
    N1 = len(all_dists)
    N2 = len(all_dists[0])
    for dists in all_dists:
        assert len(dists) == N2

    D = np.array(all_dists)
    print('loaded dists',D.shape)
    return D

def get_file1_tcr_efforts( repfile1, repfile2, verbose=True ):
    ''' Return the row sums of the effort matrix, ie the per-tcr efforts for each tcr in repfile1

    to get tcrdist-scaled output the row sums are multiplied by the size of the repertoire, un-doing the
    (uniform) weighting in the transport calc
    '''
    global lambd
    global Dmax
    D = get_raw_distance_matrix( repfile1, repfile2 )
    N1,N2 = D.shape

    D /= Dmax
    wts1 = np.ones((N1,)) / N1
    wts2 = np.ones((N2,)) / N2

    if verbose:
        print('run sinkhorn for mat',N1,N2)
        sys.stdout.flush()

    mat = ot.sinkhorn(wts1, wts2, D, lambd)

    # total effort matrix
    hadamard = np.multiply(D, mat)


    ## take row sums of the effort matrix, multiply by Dmax to get tcrdists and by N1 to remove wts
    repfile1_efforts = Dmax * N1 * hadamard.sum(axis=1)

    assert repfile1_efforts.shape[0] == N1

    return repfile1_efforts



## loop over the foreground repertoires (but actually we only do the first one, see early exit below)
## could run the full loop...


nbhd_result = defaultdict(dict)
for repfile1 in fg_repfiles:
    print(repfile1)
    tcrs = [x[:-1] for x in open(repfile1,'r')]
    N1 = len(tcrs)

    ## compute intra-repertoire distance matrix to find TCR neighborhoods
    D_11 = get_raw_distance_matrix( repfile1, repfile1 )
    subject = ntpath.basename(repfile1)
    np.savetxt(
        "/home/bolson2/sync/per_tcr/dist_matrices/" + subject + ".csv", 
        D_11,
        delimiter=",",
        fmt='%i'
    )


    ## compute per-tcr efforts against each of the other repertoires:
    fg_efforts = [ get_file1_tcr_efforts(repfile1,x) for x in fg_repfiles if x != repfile1 ]
    bg_efforts = [ get_file1_tcr_efforts(repfile1,x) for x in bg_repfiles if x != repfile1 ]


    T_fg = np.array( fg_efforts ).transpose()
    T_bg = np.array( bg_efforts ).transpose()

    assert T_fg.shape == (N1,len(fg_repfiles)-1)
    assert T_bg.shape == (N1,len(bg_repfiles))

    T_fg_means = np.mean(T_fg,axis=1)
    T_fg_stddevs = np.std(T_fg,axis=1)
    T_bg_means = np.mean(T_bg,axis=1)
    T_bg_stddevs = np.std(T_bg,axis=1)

    result = {"foreground": {"means": T_fg_means.tolist(), "std_devs": T_fg_stddevs.tolist()}, "background": {"means": T_bg_means.tolist(), "std_devs": T_bg_stddevs.tolist()}}
    with open("/home/bolson2/sync/per_tcr/empirical_fg_bg_stats.json", "w") as fp:
        json.dump(result, fp)

    nbhd_means = defaultdict(dict)
    nbhd_sds = defaultdict(list)
    for ii in range(N1):
        nbhd_means_by_cutoff = defaultdict(dict)
        for nbrcutoff in nbrcutoffs:
            ii_nbrhood_mask = ( D_11[ii,:] < nbrcutoff )

            # the efforts relative to the BG repertoires
            ii_bg_efforts = T_bg[ii,:]
            # the efforts rel to BG, summed over the nbrs (includes ii, might only be ii)
            ii_bg_nbrhood_efforts = np.sum( T_bg[ii_nbrhood_mask,:], axis=0 )

            # the efforts relative to the FG repertoires
            ii_fg_efforts = T_fg[ii,:] # the efforts rel to FG, summed over the nbrs (includes ii, might only be ii)
            ii_fg_nbrhood_efforts = np.sum( T_fg[ii_nbrhood_mask,:], axis=0 )

            T_fg_nbhd = np.array( ii_fg_nbrhood_efforts ).transpose()
            T_bg_nbhd = np.array( ii_bg_nbrhood_efforts ).transpose()

            
            neighbor_count = sum(ii_nbrhood_mask).item()
            cutoff_means = defaultdict(dict)
            cutoff_means[neighbor_count] = {
                "foreground": np.mean(T_fg_nbhd,axis=0), 
                "background": np.mean(T_bg_nbhd,axis=0),
            }

            nbhd_means_by_cutoff[nbrcutoff] = cutoff_means

            cutoff_sds = {
                "foreground": np.std(T_fg_nbhd,axis=0),
                "background": np.std(T_bg_nbhd,axis=0)
            }

            # compare fg, bg distributions of efforts
            mwu_stat, mwu_pval = mannwhitneyu( ii_bg_efforts, ii_fg_efforts )
            t_stat, t_pval     = ttest_ind( ii_bg_efforts, ii_fg_efforts )

            _, mwu_pval_nbrhood = mannwhitneyu( ii_bg_nbrhood_efforts, ii_fg_nbrhood_efforts )
            _, t_pval_nbrhood     = ttest_ind( ii_bg_nbrhood_efforts, ii_fg_nbrhood_efforts )

            # crude multiple-testing correction:
            mwu_pval *= N1
            t_pval *= N1
            mwu_pval_nbrhood *= N1
            t_pval_nbrhood *= N1

            print('fg_trans: {:6.2f} {:6.2f} bg: {:6.2f} {:6.2f} t_stat: {:6.2f} t_pval: {:8.1e} {:8.1e} mwu_pval: {:8.1e} {:8.1e} nbrs: {:3d} {}'\
                  .format( T_fg_means[ii], T_fg_stddevs[ii], T_bg_means[ii], T_bg_stddevs[ii],
                           t_stat, t_pval, t_pval_nbrhood, mwu_pval, mwu_pval_nbrhood, np.sum(ii_nbrhood_mask), tcrs[ii]))
        nbhd_means[ii] = defaultdict(dict)
        nbhd_means[ii][tcrs[ii]] = nbhd_means_by_cutoff

    nbhd_result[subject] = {"means": nbhd_means} #, "sds": nbhd_sds }
    with open("/home/bolson2/sync/per_tcr/empirical_fg_bg_nbhd_stats.json", "w") as fp:
        json.dump(nbhd_result, fp)

