## this script bundles various exploratory yfv analysis steps
## sorry for the confusing mash-up
##
## raw data comes from https://github.com/mptouzel/pogorelyy_et_al_2018.git
## which goes with the publication: https://www.pnas.org/content/115/50/12704
##
from glob import glob
import numpy as np
import matplotlib.pylab as plt
import ot
from glob import glob
from os import popen, remove
from os.path import exists
import sys
import random

basedir = '/home/pbradley/csdat/yfv/pogorelyy_et_al_2018/Yellow_fever/'


if 1: # make tree from sinkhorns

    logfile = 'results/yfv/all_sinkhorns.txt'

    ## load the distances from the logfile for step 1
    all_ots = {}
    for line in open(logfile,'r'):
        l = line.split()
        assert ':sinkhorn_dist:' in l[0]
        ll = l[0].split('_vs_')
        m1 = ll[0][: ll[0].index('.txt') ]
        m2 = ll[1][: ll[1].index('.txt') ]
        dist = float(l[1])
        all_ots[(m1,m2)] = dist
        all_ots[(m2,m1)] = dist

    nrows=1
    ncols=2
    plt.figure(1,figsize=(6*ncols,6*nrows))

    from scipy.cluster import hierarchy
    from scipy.spatial import distance

    plotno=0
    if 1:
        plotno += 1
        ax = plt.subplot(nrows,ncols,plotno)

        mice = sorted( set( [x[0] for x in all_ots ] ) )

        N = len(mice)

        D = np.zeros( (N,N) )
        for i,m1 in enumerate(mice):
            for j,m2 in enumerate(mice):
                if i<=j:continue
                dist = all_ots[(m1,m2)]
                D[i,j] = dist
                D[j,i] = dist

        Z = hierarchy.linkage( distance.squareform(D,force='tovector'), method='average' )

        hierarchy.dendrogram( Z, ax=ax, orientation='right', labels = mice )
        plt.title('OT dendrogram')


        leaves = list( hierarchy.leaves_list( Z ) )
        leaves.reverse() ## imshow order swapped
        D2 = np.zeros((N,N))
        all_vals = []
        for i,li in enumerate(leaves):
            for j,lj in enumerate(leaves):
                D2[i,j] = D[li,lj]
                all_vals.append(D2[i,j])
        all_vals.sort()
        vmin = 0.
        vmax = all_vals[ int(0.9*len(all_vals)) ]
        print( 'vmax=',vmax )

        plotno += 1
        plt.subplot(nrows,ncols,plotno)
        plt.imshow(D2,interpolation='nearest',vmin=0.,vmax=vmax)
        plt.yticks( range(N), [ mice[x] for x in leaves], fontsize=7 )
        plt.xticks( range(N), [ mice[x] for x in leaves], fontsize=7, rotation = 'vertical' )
        plt.title('OT matrix')
        plt.colorbar()



#     plt.suptitle("""Individual IELrep mouse repertoire comparisons by OT
# A is alpha chain repertoire OT
# B is beta chain repertoire OT
# AB is alpha OT + beta OT""")
    plt.subplots_adjust(top = 0.93, bottom = 0.06, hspace = 0.15 )
    tree_pngfile = 'results/yfv/OT_tree_and_matrix.png'
    print('making:',tree_pngfile)
    plt.savefig(tree_pngfile,dpi=150)

    exit()


if 0: # setup for running on cluster to compute the sinkhorns: make a long list of commands in a text file
    #
    import random
    from os.path import exists
    from sys import exit

    py_exe = '/home/pbradley/anaconda2/envs/pottest/bin/python'
    script = '/home/pbradley/tcr_scripts/find_lonely_yfv_tcrs.py'

    # this was run where the tcrs files existed, won't generalize to new transport repo but FWIW
    files = glob('{}[PQS][12]*F[12]*tcrs'.format(basedir))

    outdir = '{}run1_output/'.format(basedir)
    assert exists(outdir)

    cmds = []
    for f1 in files:
        for f2 in files:
            if f1 == f2:
                continue
            outfile = '{}{}_vs_{}.txt'.format(outdir,f1.split('/')[-1],f2.split('/')[-1])
            cmd = '{} {} {} {} > {}'.format(py_exe, script, f1, f2, outfile )
            cmds.append( cmd )

    random.shuffle(cmds)
    outfile = '{}/commands.txt'.format(outdir)
    print( 'making',outfile)
    out = open(outfile,'w')
    out.write('\n'.join(cmds)+'\n')
    out.close()
    exit()

if 0: # process the raw yfv data from pogo github to make .tcrs files
    # https://github.com/mptouzel/pogorelyy_et_al_2018.git
    #

    #['Clone ID', 'Clone count', 'Clone fraction', 'Clonal sequence(s)', 'Clonal sequence quality(s)', 'All V hits', 'All D hits', 'All J hits', 'All C hits', 'All V alignments', 'All D alignments', 'All J alignments', 'All C alignments', 'N. Seq. FR1', 'Min. qual. FR1', 'N. Seq. CDR1', 'Min. qual. CDR1', 'N. Seq. FR2', 'Min. qual. FR2', 'N. Seq. CDR2', 'Min. qual. CDR2', 'N. Seq. FR3', 'Min. qual. FR3', 'N. Seq. CDR3', 'Min. qual. CDR3', 'N. Seq. FR4', 'Min. qual. FR4', 'AA. Seq. FR1', 'AA. Seq. CDR1', 'AA. Seq. FR2', 'AA. Seq. CDR2', 'AA. Seq. FR3', 'AA. Seq. CDR3', 'AA. Seq. FR4', 'Ref. points']

    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', \
                   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']


    # make tcrs files
    files = glob('[PQS][12]*F[12]*.txt')

    topn = 1000

    for file in files:
        print(file)
        header = []
        tcrs = []
        for line in open(file,'rU'):
            l = line[:-1].split('\t')
            if not header:
                header = l[:]
                print( header )
            else:
                cdr3 = l[ header.index('AA. Seq. CDR3') ]
                badseq = len(cdr3)<=5
                for aa in cdr3:
                    if aa not in amino_acids:
                        badseq = True
                        print('bad:',cdr3)
                        break
                if badseq:
                    continue
                v_hits = l[ header.index('All V hits') ]
                vg = v_hits.split(',')[0]
                vg = vg[:vg.index('*')]+'*01'
                tcrs.append( ( vg+','+cdr3 ) )
                if len(tcrs)>=topn:
                    break

        assert len(tcrs) == topn
        outfile = '{}.top{}.tcrs'.format(file,topn)
        print('making:',outfile)
        out = open(outfile,'w')
        out.write('\n'.join(tcrs)+'\n' )
        out.close()

    exit()


## below this is the script that was run on the cluster to compute sinkhorns and nbrhood loneliness values for
## all pairs of YFV repertoires
##

## this is the tcrdist-computing executable
exe = '/home/pbradley/gitrepos/pubtcrs/bin/tcrdists'

## this is a db-directory needed for the tcrdists calc, for mouse tcrs
#db = '/loc/no-backup/pbradley/share/pot_data/fake_pubtcrs_db_mouse'
db = '/home/pbradley/gitrepos/pubtcrs/db' # human

## these are the two tcrs files
#file1 = basedir+'P1_0_F1_.txt.top1000.tcrs'
#file2 = basedir+'P1_0_F2_.txt.top1000.tcrs'

file1 = sys.argv[1]
file2 = sys.argv[2]

assert exists(file1) and exists(file2)

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
      .format(dist[0]*Dmax,lambd,Dmax,N1,N2,file1,file2))
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


