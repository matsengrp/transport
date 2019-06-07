## This script computes OT distances between the individual mouse repertoires
## in the IEL dataset and makes a dendrogram from the distance matrix
##
## It is run in two steps, since the distance computation is slow (~20 minutes on rhino1)
## and the figure-making necessitates tweaking
##
## in step1, the distances are computed and written out to a logfile
## in step2, the logfile is parsed and a tree is written out
##

import numpy as np
import matplotlib.pyplot as plt
import ot
#import ot.plot
from glob import glob
from os import popen, remove, system
from os.path import exists
import sys
import random

## this is the tcrdist-computing executable
exe = '/home/pbradley/gitrepos/pubtcrs/bin/tcrdists'

## this is a db-directory needed for the tcrdists calc, for mouse tcrs
db = '/loc/no-backup/pbradley/share/pot_data/fake_pubtcrs_db_mouse'

## this is the directory where the individual mouse tcrs files are
tcrs_dir = '/loc/no-backup/pbradley/share/pot_data/iels_tcrs_by_mouse/'

#lambd = .025
lambd = .01

Dmax = 200 # constant across comparisons

# run this script once with step=1 to compute the distances and save the output in a logfile
# then run it again with step=2 to make the tree
#
step = 2
step1_logfile = 'ot_dists_iels.txt' # wherever you put the output from step 1

# this is the name of the pngfile to be generated:
tree_pngfile = 'IELrep_OT_mouse_trees.png'

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

if step == 1:
    # compute the distances for use in step 2
    for ab in 'AB':
        tcr_files = sorted( glob('{}*_{}.tcrs'.format(tcrs_dir,ab)) )

        for f1 in tcr_files:
            for f2 in tcr_files:
                if f2<f1:
                    continue
                D = get_raw_distance_matrix( f1, f2 )
                D /= Dmax
                N1,N2 = D.shape

                wts1 = np.ones((N1,)) / N1
                wts2 = np.ones((N2,)) / N2

                print('run sinkhorn2 for dist',N1,N2)
                dist = ot.sinkhorn2(wts1, wts2, D, lambd )

                assert len(dist) == 1
                print('sinkhorn_dist: {} lambda: {} Dmax: {} Ns: {} {} files: {} {}'\
                      .format(dist[0]*Dmax,lambd,Dmax,N1,N2,f1,f2))
                sys.stdout.flush()

else:
    assert step == 2

    ## load the distances from the logfile for step 1
    all_ots = { 'A': {}, 'B': {}, 'AB': {} }
    for line in open(step1_logfile,'r'):
        l = line.split()
        if l and l[0] == 'sinkhorn_dist:':
            dist = float(l[1])
            f1,f2 = l[-2:]
            ab = f1[-6] # A or B
            m1 = f1.split('/')[-1][:-7] # eg: "CD4_17"
            m2 = f2.split('/')[-1][:-7]
            all_ots[ab][(m1,m2)] = dist
            all_ots[ab][(m2,m1)] = dist

    # define alpha+beta chain distances as just the sum of the alpha distance and the beta distance
    for m1_m2 in all_ots['A']:
        all_ots['AB'][m1_m2] = all_ots['A'][m1_m2] + all_ots['B'][m1_m2]

    nrows = 3
    ncols = 1

    plt.figure(1,figsize=(8*ncols,8*nrows))

    from scipy.cluster import hierarchy
    from scipy.spatial import distance

    plotno=0
    for ab in ['A','B','AB']: # alpha, beta, and alpha+beta OT repertoire distances
        plotno += 1
        ax = plt.subplot(nrows,ncols,plotno)

        mice = sorted( set( [x[0] for x in all_ots[ab] ] ) )

        N = len(mice)

        D = np.zeros( (N,N) )
        for i,m1 in enumerate(mice):
            for j,m2 in enumerate(mice):
                if i<=j:continue
                dist = all_ots[ab][(m1,m2)]
                D[i,j] = dist
                D[j,i] = dist

        Z = hierarchy.linkage( distance.squareform(D,force='tovector'), method='average' )

        hierarchy.dendrogram( Z, ax=ax, orientation='right', labels = mice )
        plt.title('{} OT dendrogram'.format(ab))

    plt.suptitle("""Individual IELrep mouse repertoire comparisons by OT
A is alpha chain repertoire OT
B is beta chain repertoire OT
AB is alpha OT + beta OT""")
    plt.subplots_adjust(top = 0.93, bottom = 0.03, hspace = 0.15 )
    print('making:',tree_pngfile)
    plt.savefig(tree_pngfile,dpi=150)



# if 0: # copy the tcrs files to shared location and rename
#     newdir = '/loc/no-backup/pbradley/share/pot_data/iels_tcrs_by_mouse/'
#     files = glob('/home/pbradley/csdat/tcr-dist/iels/new_191218/condIELrep_*probfilt.tcrs')
#     for file in files:
#         pos1 = file.index('_cell')
#         pos2 = file.index('_sub')
#         pos3 = file.index('_chain')
#         ct  = file[ pos1+5:pos2 ]
#         sub = file[ pos2+4:pos3 ]
#         ab  = file[ pos3+6 ]
#         mouse = int(sub.split('.')[1] )
#         newfile = '{}{}_{}_{}.tcrs'.format(newdir,ct,mouse,ab)
#         system('cp {} {}'.format(file,newfile))

#     exit()

