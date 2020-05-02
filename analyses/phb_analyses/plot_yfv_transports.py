## this script reads pre-computed loneliness values that were generated by find_lonely_yfv_tcrs.py
##

from __future__ import print_function
from glob import glob
import matplotlib
from os import popen
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

rundir = 'data/yfv_raw/' #'/loc/no-backup/pbradley/share/pot_data/yfv/'

# files look like:
# Q1_0_F1_.txt.top1000.tcrs_vs_Q1_0_F2_.txt.top1000.tcrs.txt
files = glob('{}*_vs_*.txt'.format(rundir))


all_vals = {0:[], 1:[]}

for file in files:
    ll = file.split('/')[-1].split('_vs_')
    s1 = ll[0][: ll[0].index('_.txt') ].split('_') # eg ['P2','0','F1']
    s2 = ll[1][: ll[1].index('_.txt') ].split('_')
    assert len(s1) == 3 and len(s2) == 3
    if s1[0] != s2[0]:
        # different subjects
        continue
    elif s1[1] == s2[1]:
        selfcomp = 1
    elif s1[1] == '15' and s2[1] == '0':
        selfcomp = 0
    else:
        continue

    print( selfcomp, s1, s2, file )
    # read the loneliness values for the "s1" tcrs versus repertoire "s2"
    for line in popen('grep ^f1_ind ' + file ):
        l = line.split()
        all_vals[selfcomp].append( ( float(l[3]), float(l[5]) ) )

plt.figure(1,figsize=(12,12))

## scatter plots
for ii in range(2):
    plt.subplot(2,2,ii+1)
    xmx = max( x[0] for y in all_vals.values() for x in y )
    ymx = max( x[1] for y in all_vals.values() for x in y )
    plt.scatter( [x[0] for x in all_vals[ii] ], [x[1] for x in all_vals[ii]], alpha=0.33 )
    plt.xlim((-5,xmx+10))
    plt.ylim((-25,ymx+10))
    plt.xlabel('TCR self-transport')
    plt.ylabel('TCR neighborhood-transport')
    if ii==0:
        plt.title('timepoint 15 versus timepoint 0 (total TCRs= {})'.format(len(all_vals[ii])))
    else:
        plt.title('biological replicate 1 vs replicate 2 (total TCRs= {})'.format(len(all_vals[ii])))


## KDE plots
for jj in range(2):
    plt.subplot(2,2,jj+3)
    for ii in range(2):
        maxval = max( x[jj] for x in all_vals[ii] )
        xs = np.linspace(0,maxval,100)
        dens = gaussian_kde( [x[jj] for x in all_vals[ii]] )
        ys = dens(xs)
        plt.plot( xs, dens(xs), label = '15-vs-0 rep1-vs-rep2'.split()[ii])
    plt.legend()
    if jj==0:
        plt.title('Gaussian KDE: TCR self-transport')
    else:
        plt.title('Gaussian KDE: TCR neighborhood-transport')


plt.suptitle("""Comparison of loneliness measures: "self-transport" (row/col sum of Hadamard matrix) versus
"neighborhood-transport" (sum of self-transports for all TCRs within TCRdist of 48.5)""")

plt.subplots_adjust(top=0.9,bottom=0.05)
plt.savefig('results/yfv/yfv_transports.png')

