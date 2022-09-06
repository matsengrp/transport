## This is the script that we used to setup and run the ALICE spike-in experiment
##
## To get this to run on your system, you will need to
## 1. download the ALICE github repository from https://github.com/pogorely/ALICE
## 2. copy the run_alice_100M.R script to the ALICE/ repo folder
## 3. modify that run_alice_100M.R script to fixup the path for your system
##     (search for pbradley there)
## 4. modify this script (wherever it says 'pbradley') to change the alice_dir
##    variable to point to an output folder for the ALICE results, and the
##    scriptlines variable to point to the location of your ALICE repo and
##    remove or modify the "ml" command in there to whatever has to happen in your
##    system to make the Rscript command available (we use ml aka module load)
## 5. Run the script once to generate the sbatch input file
## 6. Run the sbatch commands in "cmds_file" below to launch the ALICE jobs on your
##    computing cluster.
## 7. Re-run the script with the "if 0" changed to "if 1" to analyze the results.
## 8. email pbradley@fredhutch.org if you run into any trouble or have questions!
##
import numpy as np
import ot
from glob import glob
from os import popen, remove, system
from os.path import exists
import sys
import random
import pandas as pd

naive_tcrs_file = 'data/flu_spike_ins/naive_tcrs.tsv'
flu_tcrs_file = 'data/flu_spike_ins/flu_tcrs.tsv'

random_state = 10

num_naive_list = [1000]
num_flu_list = [5, 10, 20, 40, 80, 160, 320]
num_repeats = 10

cdr3_col = 'CDR3.amino.acid.sequence'
min_nbrs = 3

# change this, it's the place where the ALICE results get written:
alice_dir = '/home/pbradley/csdat/transport/flu_spike_ins/alice/'

######################################################################################88
######################################################################################88
######################################################################################88

if 0: # read the results
    # change this to "if 1" after the ALICE calculations have been run
    results_file = 'results/flu_spike_ins/flu_auc_values_alice_run2.tsv'
    #results_file = 'results/flu_spike_ins/flu_auc_values_alice.tsv'
    runtag = 'run2' # 100M TCRs
    #runtag = 'run1'
    naive_tcrs_df = pd.read_table(naive_tcrs_file)
    flu_tcrs_df = pd.read_table(flu_tcrs_file)

    dfl = []
    for num_naive in num_naive_list:
        for num_flu in num_flu_list:
            for r in range(num_repeats):
                naive_sample = naive_tcrs_df.sample(
                    2*num_naive-num_flu, random_state=random_state+r)
                fg_df = pd.concat([
                    flu_tcrs_df.sample(num_flu, random_state=random_state+r),
                    naive_sample.iloc[num_naive:],
                ])
                cols = 'va ja cdr3a vb jb cdr3b'.split()
                fg_tcrs = list(fg_df[cols].itertuples(index=False, name=None))
                files = glob(f'{alice_dir}N{num_naive}_F{num_flu}_R{r}_*{runtag}.csv')
                if files:
                    df = pd.concat([pd.read_csv(x) for x in files])
                    num_hits = df.shape[0]
                    #print(df.iloc[0])
                    tcrs = list(df[cols].itertuples(index=False, name=None))
                    inds = [fg_tcrs.index(x) for x in tcrs]
                    flu_count = sum(x<1000 for x in inds)
                    naive_count = sum(x>=1000 for x in inds)
                else:
                    flu_count, naive_count = 0, 0
                flu_total = num_flu
                naive_total = num_naive - num_flu
                assert naive_count == 0
                assert flu_count <= flu_total
                #auc = ( (naive_count/naive_total)*(0.5*(0+flu_count/flu_total)) +
                #        (...  ### more complicated if naive_count>0
                auc = 0.5 * (flu_count/flu_total + 1.0)
                #print(num_naive, num_flu, r, auc, flu_count, naive_count)

                dfl.append(dict(
                    num_naive=num_naive,
                    num_flu=num_flu,
                    flu_count=flu_count,
                    naive_count=naive_count,
                    auc_alice=auc,
                    repeat=r,
                ))


    pd.DataFrame(dfl).to_csv(results_file, sep='\t', index=False)
    print('made:', results_file)
    exit()


## change this to point to your ALICE repo
scriptlines = """#!/bin/bash

ml fhR
cd /home/pbradley/gitrepos/ALICE/

Rscript run_alice_100M.R $1 $2 $3
"""

# $1 is the tsv file
# $2 is the output folder for the .rda files
# $3 is the output csv file

min_nbrs = 3

runtag = 'run2'
#runtag = 'run1'
cmds_file = f'{alice_dir}{runtag}_commands.txt'
assert not exists(cmds_file)
out = open(cmds_file, 'w')

script = f'{alice_dir}{runtag}_script.sh'
with open(script, 'w') as f:
    f.write(scriptlines)
system('chmod a+x '+script)


naive_tcrs_df = pd.read_table(naive_tcrs_file)
flu_tcrs_df = pd.read_table(flu_tcrs_file)


def count_mismatches(a,b):
    return sum(x!=y for x,y in zip(a,b))

dfl = []
counter=0
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

            # look for potential alice nbrhoods
            # ALICE only works for beta chains
            for ab in ['b']: #'ab':
                df = fg_df.copy()
                df['Read.count'] = 10
                df['bestVGene'] = [x[:x.index('*')] for x in df['v'+ab]]
                df['bestJGene'] = [x[:x.index('*')] for x in df['j'+ab]]
                df[cdr3_col] = df['cdr3'+ab]
                df['cdr3len'] = df[cdr3_col].str.len()
                df['vjl'] = df.bestVGene+'_'+df.bestJGene+'_'+df.cdr3len.astype(str)
                vj_combos = set()
                for vjl, count in df.vjl.value_counts().items():
                    has_nbrs = False
                    if count < min_nbrs+1:
                        break
                    cdr3s = list(df[df.vjl==vjl][cdr3_col])
                    for i, a in enumerate(cdr3s):
                        nbrs=0
                        for j, b in enumerate(cdr3s):
                            if i!=j:
                                nbrs += count_mismatches(a,b)<=1
                        if nbrs >= min_nbrs:
                            #print('nbrs:', a, nbrs, vjl, tsvfile)
                            has_nbrs = True
                    if has_nbrs:
                        vj_combos.add(tuple(vjl.split('_')[:2]))
                for v, j in vj_combos:
                    # make tsv file with just the vj tcrs
                    # add new command line
                    df_vj = df[(df.bestVGene==v)&(df.bestJGene==j)]
                    vj_tag = (v.replace('-','_').replace('/','_')+'_'+
                              j.replace('-','_').replace('/','_'))
                    out_tsvfile = (f'{alice_dir}N{num_naive}_F{num_flu}_R{r}_'
                                   f'{ab}_{vj_tag}.tsv')
                    df_vj.to_csv(out_tsvfile, sep='\t', index=False)
                    counter += 1

                    outprefix = out_tsvfile[:-4]+'_'+runtag
                    cmd = (f'sbatch -t 1-0 -e {outprefix}.err -o {outprefix}.log '
                           f' {script} {out_tsvfile} folder{counter} {outprefix}.csv')
                    out.write(cmd+'\n')



out.close()
print('made:', cmds_file)

