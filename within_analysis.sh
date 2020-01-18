#!/bin/sh
METHOD=$1

rm -rf ~/sync/$METHOD/*

python within_gene.py CD4_17 $METHOD
Rscript --vanilla plot_within_results.R CD4 $METHOD

python within_gene.py CD8_10 $METHOD
Rscript --vanilla plot_within_results.R CD8 $METHOD

python within_gene.py DN_15 $METHOD
Rscript --vanilla plot_within_results.R DN $METHOD

python permutation.py $METHOD

cp -r $METHOD/* ~/sync/$METHOD/

