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
exe = '/loc/no-backup/pbradley/share/pot_data/tcrdists'

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
step = 1
step1_logfile = 'results/iels/ot_dists_iels.txt' # wherever you put the output from step 1

# this is the name of the pngfile to be generated:
tree_pngfile = 'results/iels/IELrep_OT_mouse_trees.png'
single_tree_pngfile = 'results/iels/IELrep_OT_mouse_tree_beta.png'
#single_tree_pngfile = '/home/pbradley/csdat/tmp.png'

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

if 0:
    # look for overlap between CD8 and DN
    ab = 'B'
    reps = 'CD4 CD8 DN'.split()
    num_mice = 23
    for r0 in reps:
        for r1 in reps:
            if r0 <= r1: continue
            for ii in range(1,num_mice+1):
                f0 = f'{tcrs_dir}{r0}_{ii}_{ab}.tcrs'
                f1 = f'{tcrs_dir}{r1}_{ii}_{ab}.tcrs'
                assert exists(f0) and exists(f1)
                tcrs0 = set([x.strip() for x in open(f0,'r').readlines()])
                tcrs1 = set([x.strip() for x in open(f1,'r').readlines()])
                overlap = len(tcrs0&tcrs1)
                jac = overlap/len(tcrs0|tcrs1)
                print(f'jaccard: {jac:.6f} {overlap:2d} {r0:3s} {r1:3s} mouse: {ii:2d}')

    exit()


if step == 1:
    import time
    # compute the distances for use in step 2
    tcrdist_time = 0
    sinkhorn_time = 0
    for ab in 'AB':
        tcr_files = sorted( glob('{}*_{}.tcrs'.format(tcrs_dir,ab)) )

        for f1 in tcr_files:
            for f2 in tcr_files:
                if f2<f1:
                    continue
                t0 = time.time()
                D = get_raw_distance_matrix( f1, f2 )
                tcrdist_time += time.time() - t0
                D /= Dmax
                N1,N2 = D.shape

                wts1 = np.ones((N1,)) / N1
                wts2 = np.ones((N2,)) / N2

                print('run sinkhorn2 for dist',N1,N2)
                t0 = time.time()
                dist = ot.sinkhorn2(wts1, wts2, D, lambd, verbose=True)
                if type(dist) is list or type(dist) is np.ndarray:
                    dist = dist[0]
                print(dist)
                sinkhorn_time += time.time() - t0
                print('sinkhorn_dist: {} lambda: {} Dmax: {} Ns: {} {} files: {} {}'\
                      .format(dist*Dmax,lambd,Dmax,N1,N2,f1,f2))
                print(f'times: tcrdist: {tcrdist_time} sinkhorn: {sinkhorn_time}')
                sys.stdout.flush()

else:
    assert step == 2

    ## load the distances from the logfile for step 1
    all_ots = { 'A': {}, 'B': {}, 'AB': {} }
    all_tcr_counts = { 'A': {}, 'B': {} }
    for line in open(step1_logfile,'r'):
        l = line.split()
        if l and l[0] == 'sinkhorn_dist:':
            assert l[6] == 'Ns:'
            dist = float(l[1])
            f1,f2 = l[-2:]
            N1,N2 = int(l[7]), int(l[8])
            ab = f1[-6] # A or B
            m1 = f1.split('/')[-1][:-7] # eg: "CD4_17"
            m2 = f2.split('/')[-1][:-7]
            all_ots[ab][(m1,m2)] = dist
            all_ots[ab][(m2,m1)] = dist
            all_tcr_counts[ab][m1] = N1
            all_tcr_counts[ab][m2] = N2

    assert sorted(all_tcr_counts['A'].keys()) == sorted(all_tcr_counts['B'].keys())
    all_tcr_counts['AB'] = {
        m:all_tcr_counts['A'][m]+all_tcr_counts['B'][m]
        for m in all_tcr_counts['A']
    }
    # define alpha+beta chain distances as just the sum of the alpha distance and the beta distance
    for m1_m2 in all_ots['A']:
        all_ots['AB'][m1_m2] = all_ots['A'][m1_m2] + all_ots['B'][m1_m2]

    nrows = 3
    ncols = 1

    plt.figure(1,figsize=(8*ncols,8*nrows))

    from scipy.cluster import hierarchy
    from scipy.spatial import distance

    plotno=0
    save_info = {}
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

        Z = hierarchy.linkage( distance.squareform(D,force='tovector'), method='ward' )

        labels = [f'{m} (#={all_tcr_counts[ab][m]})' for m in mice]

        hierarchy.dendrogram( Z, ax=ax, orientation='right', labels = labels )
        plt.title('{} OT dendrogram'.format(ab))

        save_info[ab] = [Z, mice]

    plt.suptitle("""Individual IELrep mouse repertoire comparisons by OT
A is alpha chain repertoire OT
B is beta chain repertoire OT
AB is alpha OT + beta OT""")
    plt.subplots_adjust(top = 0.93, bottom = 0.03, hspace = 0.15 )
    print('making:',tree_pngfile)
    plt.savefig(tree_pngfile,dpi=150)


    # just make a single tree, demo for figure
    plt.figure(figsize=(8,8))
    rep_color = {'CD4':'C0', 'CD8':'C1', 'DN':'C2'}
    Z, mice = save_info['B']
    hierarchy.dendrogram(Z, ax=plt.gca(), orientation='right', labels = mice,
                         link_color_func=lambda x:'gray')##AAAAAA')
    #for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), colors):
    for ticklabel in plt.gca().get_yticklabels():
        color = rep_color[ticklabel.get_text().split('_')[0]]
        print(ticklabel)
        #print(dir(ticklabel))
        #exit()
        ticklabel.set_color(color)


    # locs, labels = plt.yticks()
    # print(locs)
    # print(labels)
    # leaves = hierarchy.leaves_list(Z)
    # labels = [mice[x] for x in leaves]
    # colors = [rep_color[m.split('_')[0]] for m in labels]
    # plt.yticks(range(len(labels)), labels)
    # for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), colors):
    #     ticklabel.set_color(tickcolor)
    print('making:',single_tree_pngfile)
    plt.tight_layout()
    plt.savefig(single_tree_pngfile, dpi=150)




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

