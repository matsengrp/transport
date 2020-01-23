## this script reads two tcr repertoires, takes random subsets of each,
## computes tcrdist distances on the fly and the ot transport matrix
##
## then it computes for each tcr in repertoire1 the total transport in a
## nbrhood around that tcr
##

import numpy as np
#import matplotlib.pylab as pl
import ot
#import ot.plot
from glob import glob
from os import popen, remove
from os.path import exists
import sys
import random

## this is the tcrdist-computing executable
exe = '/home/pbradley/gitrepos/pubtcrs/bin/tcrdists'
#exe = 'bin/tcrdists'

## this is a db-directory needed for the tcrdists calc, for mouse tcrs
db = '/loc/no-backup/pbradley/share/pot_data/fake_pubtcrs_db_mouse'
#db = 'data/databases/fake_pubtcrs_db_mouse'

## these are the two tcrs files
file1 = '/loc/no-backup/pbradley/share/pot_data/ielrep_beta_DN_tcrs.txt'
file2 = '/loc/no-backup/pbradley/share/pot_data/ielrep_beta_CD4_tcrs.txt'
#file1 = 'data/iel_data/ielrep_beta_DN_tcrs.txt'
#file2 = 'data/iel_data/ielrep_beta_CD4_tcrs.txt'

#lambd = .025
lambd = .01

nbrcutoff = 48.5 # distance at which two single-chain tcrs are considered nbrs

Dmax = 200 # constant across comparisons

## HACK: we make temporary files with random subsets of the full datasets, then call the cmdline distance calc
tmpfile_prefix = 'tmptcrs.{}'.format(random.random())


# choose this many of each dataset
N1, N2 = 1000, 1001
#N1, N2 = 10000, 8000

cleanup_tmpfiles = True

def get_raw_distance_matrix( f1, f2 ):
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

def make_random_subset_file_and_return_lines( oldfile, N, newfile ):
    lines = open(oldfile,'r').readlines()
    random.shuffle( lines )
    lines = lines[:N]
    out = open(newfile,'w')
    out.writelines(lines)
    out.close()
    return lines

# the tmpfiles
f1 = '{}_f1_tcrs_subset.txt'.format(tmpfile_prefix)
f2 = '{}_f2_tcrs_subset.txt'.format(tmpfile_prefix)

lines1 = make_random_subset_file_and_return_lines( file1, N1, f1 )
lines2 = make_random_subset_file_and_return_lines( file2, N2, f2 )

N1 = len(lines1) # in case there are actually fewer than we requested
N2 = len(lines2) # in case there are actually fewer than we requested

assert exists(f1) and exists(f2)

D_11 = get_raw_distance_matrix( f1, f1 )
D_12 = get_raw_distance_matrix( f1, f2 )
D_22 = get_raw_distance_matrix( f2, f2 )

assert D_11.shape == (N1,N1)
assert D_12.shape == (N1,N2)
assert D_22.shape == (N2,N2)

D = np.array( D_12 ) # make copy, then normalize
D /= Dmax
wts1 = np.ones((N1,)) / N1
wts2 = np.ones((N2,)) / N2

print('run sinkhorn2 for dist',N1,N2)
dist = ot.sinkhorn2(wts1, wts2, D, lambd )

# is this silly to run again?
print('run sinkhorn for mat',N1,N2)
mat = ot.sinkhorn(wts1, wts2, D, lambd)

# total effort matrix or something
hadamard = np.multiply(D, mat)


assert len(dist) == 1
print('sinkhorn_dist: {} lambda: {} Dmax: {} Ns: {} {} files: {} {}'\
      .format(dist[0]*Dmax,lambd,Dmax,N1,N2,f1,f2))
sys.stdout.flush()


## take row and column (unused) sums of the effort matrix
f1_trans = hadamard.sum(axis=1)
f2_trans = hadamard.sum(axis=0)
assert f1_trans.shape[0] == N1
assert f2_trans.shape[0] == N2


# look through the file1 tcrs
# get total transport in nbrhood
# write out

for ii in range(N1):
    nbr_trans = []
    for jj in range(N1):
        if D_11[ii,jj] <= nbrcutoff:
            nbr_trans.append( f1_trans[jj] )
    tot = sum( nbr_trans )
    print( 'f1_index: {:4d} self_transport: {:9.3f} nbrhood_transport: {:9.3f} nbrhood_size: {:2d} tcr: {}'\
           .format( ii, Dmax * f1_trans[ii] * N1, Dmax * tot * N1, len(nbr_trans), lines1[ii][:-1] ) )


if cleanup_tmpfiles:
    remove( f1 )
    remove( f2 )


