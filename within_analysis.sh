#!/bin/sh

python within_gene.py CD4_22
Rscript --vanilla plot_within_results.R CD4

python within_gene.py CD8_19
Rscript --vanilla plot_within_results.R CD8

python within_gene.py DN_2
Rscript --vanilla plot_within_results.R DN
