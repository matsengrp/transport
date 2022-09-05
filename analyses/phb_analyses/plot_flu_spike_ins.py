#
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

results_file = 'results/flu_spike_ins/flu_auc_values.tsv'
single_chain_results_file = 'results/flu_spike_ins/flu_auc_values_single_chain.tsv'
alice_results_file = 'results/flu_spike_ins/flu_auc_values_alice_run2.tsv' # 100M -- results are the same...
#pngfile = '/home/pbradley/csdat/transport/flu_spike_in_results_run2.pdf'
#alice_results_file = 'results/flu_spike_ins/flu_auc_values_alice.tsv'
pngfile = '/home/pbradley/csdat/transport/flu_spike_in_results.pdf'
supp_pngfile = '/home/pbradley/csdat/transport/flu_spike_in_results_supp.pdf'
#pngfile = '/home/pbradley/csdat/transport/tmp.png'
#pngfile = 'results/flu_spike_ins/flu_auc_boxes.png'

index_cols = 'num_naive num_flu repeat'.split()
results = pd.read_table(results_file)
sc_results = pd.read_table(single_chain_results_file)
alice_results = pd.read_table(alice_results_file).set_index(index_cols)

results = results.join(sc_results[sc_results.chain=='A'].set_index(index_cols),
                       on = index_cols, rsuffix='_alpha')
results = results.join(sc_results[sc_results.chain=='B'].set_index(index_cols),
                       on = index_cols, rsuffix='_beta')
results = results.join(alice_results, on=index_cols)
results.sort_values(['num_flu','repeat'], inplace=True)

num_flu_list = sorted(set(results.num_flu))
num_repeat = max(results.repeat)+1

assert results.shape[0] == len(num_flu_list)*num_repeat

nrows, ncols, plotno = 1,2,0
plt.figure(figsize=(ncols*4, nrows*4))
#nrows, ncols, plotno = 1,4,0
#['auc_transport', 'TotalLoneliness'],
# for col, colname in [['auc_loneliness','RelativeLoneliness'],
#                      ['auc_loneliness_alpha','RelativeLoneliness (alpha chain)'],
#                      ['auc_loneliness_beta','RelativeLoneliness (beta chain)'],
#                      ['auc_alice', 'ALICE (beta chain)']]:
for col, colname in [['auc_loneliness_beta','Transport'],
                     ['auc_alice', 'ALICE']]:
    plotno += 1
    plt.subplot(nrows, ncols, plotno)
    vals = np.array(results[col]).reshape(-1,num_repeat).T
    #print(vals)
    plt.boxplot(x=vals) # boxplot shows each column as a separate box
    if plotno==1:
        plt.ylabel(f'Recovery of spike-in TCRs (AUROC)')
    plt.xlabel('Number of spike-in TCRs (out of 1000)')
    plt.title(f'{colname}')
    #plt.title(f'Recovery of spike-in TCRs\nwhen sorting by {colname}\n({num_repeat} random repeats)')
    locs,_ = plt.xticks()
    plt.xticks(locs, [str(x) for x in num_flu_list])
    plt.ylim([0.25,1])

plt.tight_layout()
plt.savefig(pngfile)
print('made:', pngfile)


### supp fig:
nrows, ncols, plotno = 1,2,0
plt.figure(figsize=(ncols*4, nrows*4))
#nrows, ncols, plotno = 1,4,0
#['auc_transport', 'TotalLoneliness'],
# for col, colname in [['auc_loneliness','RelativeLoneliness'],
#                      ['auc_loneliness_alpha','RelativeLoneliness (alpha chain)'],
#                      ['auc_loneliness_beta','RelativeLoneliness (beta chain)'],
#                      ['auc_alice', 'ALICE (beta chain)']]:
for col, colname in [['auc_loneliness_beta','NeighborhoodLoneliness'],
                     ['auc_transport_beta', 'IndividualLoneliness']]:
    plotno += 1
    plt.subplot(nrows, ncols, plotno)
    vals = np.array(results[col]).reshape(-1,num_repeat).T
    #print(vals)
    plt.boxplot(x=vals) # boxplot shows each column as a separate box
    if plotno==1:
        plt.ylabel(f'Recovery of spike-in TCRs (AUROC)')
    plt.xlabel('Number of spike-in TCRs (out of 1000)')
    plt.title(f'{colname}')
    #plt.title(f'Recovery of spike-in TCRs\nwhen sorting by {colname}\n({num_repeat} random repeats)')
    locs,_ = plt.xticks()
    plt.xticks(locs, [str(x) for x in num_flu_list])
    plt.ylim([0.25,1])

plt.tight_layout()
plt.savefig(supp_pngfile)
print('made:', supp_pngfile)

