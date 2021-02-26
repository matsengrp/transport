Analysis scripts
================

This directory contains Python scripts which perform various analyses on the IEL and YFV datasets.
`replicates.py` is a good starting point for the IEL data as the loneliness scores computed in this script are stored and used in several other analyses, e.g. `z_scores.py`.
`yfv.py` is a good starting point for the YFV data as the results of this script are used in `query_validation_tcrs.py`.

### cluster_iels.py
This script queries every TCR in each of the DN, CD4, and CD8 repertoires against the Tremont, Revere, and Ida clusters computed in `combined_replicates.py`.
The result are csv files (found in `output/iel_clusters`) for each dataset which label each TCR as belonging to either of these three clusters, or "N/A" otherwise.


### combined_replicates.py
This script constructs the full combined DN and CD4 datasets (by concatenating all of the respective replicate datasets within each group) and computes the top lonely clusters between them.
The result are directories of cluster-related files for each cluster, and can be found in `output/hmm/cd4_dn`.


### full_dist_matrix.py
This script computes the full distance matrix between all DN TCRs to all other DN TCRs, across all of the DN subjects, and saves it to `all_subjects.csv`.
This matrix mainly used for visualizations, e.g. examining MDS plots.


### motif.py
This script operates on the loneliness scores obtained in `replicates.py`, computing the mean annulus loneliness values used to estimate the breakpoint radius of a clusters.
The result is a csv file (found at `output/csv/motif.csv`) containing mean annulus loneliness values (called `annulus_enrichment`s in the file) for each radius for each DN subject.


### query_validation_tcrs.py
This script validates the analysis performed in `yfv.py` by querying two independent datasets, one of YFV TCRs and a control dataset of CMV TCRs, against each of the computed clusters.
The result is a csv file (found at `output/hmm/yfv_hits.csv`) containing the number of hits, or TCRs residing in the cluster, for each of the top 10 clusters for each repertoire comparison.


### replicates.py
This script iterates over every DN repertoire, and computes its "foreground" loneliness scores with respect to every CD4 repertoire, as well as its "background" loneliness scores with respect to every other DN repertoire.
This script also computes z-scores between each (DN, DN) and (DN, CD4) pair which are used to benchmark the randomization z-scores obtained in `z_scores.py`.
The results can be found at `output/json/empirical_fg_bg_nbhd_stats.json` and `output/json/replicate_z_scores.json`.


### yfv.py
This script analyses the YFV data, computing the top lonely clusters for four different timepoints, 0d vs -7d, 0d vs 0d, 0d vs 15d, and 0d vs 45d, for each of the six subjects, P1, P2, S1, S2, Q1, and Q2.
The resultant directories can be found in `output/hmm`. 


### z_scores.py
This script performs a randomization test between the biggest DN TRB repertoire (subject DN_15) and all of the CD4 repertoires.
Each randomization test generates per-TCR z-scores, which are compared to the z-scores obtained in `replicates.py`. 
The main results can be found at `output/json/rand_z_scores.json`. 
