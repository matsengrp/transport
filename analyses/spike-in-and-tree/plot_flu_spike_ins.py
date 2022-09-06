## this script would be run from the transport/ directory after the transport and
## ALICE calculations had been run
##
##
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

single_chain_results_file = 'results/flu_spike_ins/flu_auc_values_single_chain.tsv'
alice_results_file = 'results/flu_spike_ins/flu_auc_values_alice_run2.tsv' # 100M
pngfile = '/home/pbradley/csdat/transport/flu_spike_in_results.pdf'
supp_pngfile = '/home/pbradley/csdat/transport/flu_spike_in_results_supp.pdf'

index_cols = 'num_naive num_flu repeat'.split()
sc_results = pd.read_table(single_chain_results_file)
alice_results = pd.read_table(alice_results_file).set_index(index_cols)

results = sc_results[sc_results.chain=='B']

results = results.join(sc_results[sc_results.chain=='A'].set_index(index_cols),
                       on=index_cols, rsuffix='_alpha')
results = results.join(alice_results, on=index_cols)
results.sort_values(['num_flu','repeat'], inplace=True)

num_flu_list = sorted(set(results.num_flu))
num_repeat = max(results.repeat)+1

assert results.shape[0] == len(num_flu_list)*num_repeat

nrows, ncols, plotno = 1,2,0
plt.figure(figsize=(ncols*4, nrows*4))
for col, colname in [['auc_loneliness','Transport'],
                     ['auc_alice', 'ALICE']]:
    plotno += 1
    plt.subplot(nrows, ncols, plotno)
    vals = np.array(results[col]).reshape(-1,num_repeat).T
    plt.boxplot(x=vals) # boxplot shows each column as a separate box
    if plotno==1:
        plt.ylabel(f'Recovery of spike-in TCRs (AUROC)')
    plt.xlabel('Number of spike-in TCRs (out of 1000)')
    plt.title(f'{colname}')
    locs,_ = plt.xticks()
    plt.xticks(locs, [str(x) for x in num_flu_list])
    plt.ylim([0.25,1])

plt.tight_layout()
plt.savefig(pngfile)
print('made:', pngfile)


### supp fig:
nrows, ncols, plotno = 1,2,0
plt.figure(figsize=(ncols*4, nrows*4))
for col, colname in [['auc_loneliness','NeighborhoodLoneliness'],
                     ['auc_transport', 'IndividualLoneliness']]:
    plotno += 1
    plt.subplot(nrows, ncols, plotno)
    vals = np.array(results[col]).reshape(-1,num_repeat).T
    plt.boxplot(x=vals) # boxplot shows each column as a separate box
    if plotno==1:
        plt.ylabel(f'Recovery of spike-in TCRs (AUROC)')
    plt.xlabel('Number of spike-in TCRs (out of 1000)')
    plt.title(f'{colname}')
    locs,_ = plt.xticks()
    plt.xticks(locs, [str(x) for x in num_flu_list])
    plt.ylim([0.25,1])

plt.tight_layout()
plt.savefig(supp_pngfile)
print('made:', supp_pngfile)

