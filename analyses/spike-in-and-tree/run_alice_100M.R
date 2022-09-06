#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args[1] is the sample tsvfile

# args[2] is the foldername

# args[3] is the output csv filename


# I'm not an R expert, but I think this next commented out part needs to happen only once,
# for example the first time this script is run, to install the packages, then it can be
# commented out
# 
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(version = "3.15")
#BiocManager::install(c("Biostrings"))
#install.packages("igraph")
#install.packages("data.table")
#install.packages("stringdist")

library(data.table)
S1d0<-fread(args[1])
S1<-list(d0=S1d0)
source("/home/pbradley/gitrepos/ALICE/ALICE.R")
S1_alice<-ALICE_pipeline(DTlist=S1,folder=args[2],cores=1,iter=50,nrec=2e6)
write.csv(S1_alice$d0, args[3])
