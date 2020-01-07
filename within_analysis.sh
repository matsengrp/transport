#!/bin/sh

rm -rf ~/sync/subsample/*

python within_gene.py CD4_15
Rscript --vanilla plot_within_results.R CD4

python within_gene.py CD8_8
Rscript --vanilla plot_within_results.R CD8

python within_gene.py DN_8
Rscript --vanilla plot_within_results.R DN
