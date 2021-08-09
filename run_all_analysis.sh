#!/usr/bin/env bash
set -eu

rm -rf tmp_output
rm -rf output

# Computes fg and bg loneliness scores, and replicate z-scores
python analyses/replicates.py  

# Computes the 3 large combined repertoires and the top lonely clusters
python analyses/combined_replicates.py  

python analyses/cluster_iels.py  # Computes cluster memberships of all DN, CD4, and CD8 repertoires based on the 3 OT clusters identified by analyses/combined_replicates.py
python analyses/z_scores.py  # Computes randomization z-scores
python analyses/motif.py  # Computes clusters for differenct radii over IEL repertoires, to be used in R/motif.R
python analyses/yfv.py  # Computes many comparisons across donors/timepoints for yfv analysis
python analyses/query_validation_tcrs.py  # Queries the yfv/cmv validation datasets against our clusters computed in analyses/yfv.py

Rscript R/load_score_datasets.R  # Loads fg_dat and bg_dat based on the output of analyses/replicates.py
Rscript R/cluster.R  # Creates plots based on the OT clusters identified in analyses
Rscript R/make_gene_plots.R  # Creates V gene frequency bars (see Fig 4 in manuscript)
Rscript R/z_score.R  # Plots replicate vs randomization z-scores
Rscript R/motif.R  # Plots mean annulus loneliness vs cluster radius (Fig 8 in manuscript)
Rscript R/yfv.R  # Analyzes and plots results of analyses/yfv.py and analyses/query_validation_tcrs.py

rm -rf tmp_output
