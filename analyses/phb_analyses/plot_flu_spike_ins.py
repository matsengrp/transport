#
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

results_file = 'results/flu_spike_ins/flu_auc_values.tsv'
#pngfile = '/home/pbradley/csdat/tmp.png'
pngfile = 'results/flu_spike_ins/flu_auc_boxes.png'

results = pd.read_table(results_file)
results.sort_values(['num_flu','repeat'], inplace=True)

num_flu_list = sorted(set(results.num_flu))
num_repeat = max(results.repeat)+1

assert results.shape[0] == len(num_flu_list)*num_repeat
plt.figure(figsize=(7.5,10))

nrows, ncols, plotno = 2,1,0
for col, colname in [['auc_loneliness','RelativeLoneliness'],
                     ['auc_transport', 'TotalLoneliness']]:
    plotno += 1
    plt.subplot(nrows, ncols, plotno)
    vals = np.array(results[col]).reshape(-1,num_repeat).T
    #print(vals)
    plt.boxplot(x=vals) # boxplot shows each column as a separate box
    plt.ylabel('AUROC')
    plt.xlabel('Number of Flu M158 TCRs spiked into repertoire of size 1000')
    plt.title(f'Recovery of spike-in TCRs when sorting by {colname}\n({num_repeat} random repeats)')
    locs,_ = plt.xticks()
    plt.xticks(locs, [str(x) for x in num_flu_list])
    plt.ylim([0.25,1])

plt.tight_layout()
plt.savefig(pngfile)

