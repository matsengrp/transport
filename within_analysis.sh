#!/bin/sh

rm -rf ~/sync/subsample/*

python within_gene.py CD4_20 ratio
Rscript --vanilla plot_within_results.R CD4 subsample

python within_gene.py CD8_15 ratio
Rscript --vanilla plot_within_results.R CD8 subsample

python within_gene.py DN_10 ratio
Rscript --vanilla plot_within_results.R DN subsample
