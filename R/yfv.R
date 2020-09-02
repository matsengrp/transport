library(data.table)
library(dplyr)

discard_extra_gene_calls <- function(genes) {
    top_gene <- genes %>% 
        sapply(strsplit, split=",") %>% 
        unlist %>% 
        first
    return(top_gene)
}

label_to_dataset <- function(label) {
    info_list <- label %>% strsplit(split="_") %>% unlist
    prefix <- info_list[1]
    timepoint <- info_list[3]
    if(timepoint == 0) {
        midfix <- "0_F2"
    } else {
        midfix <- paste(timepoint, "F1", sep="_")
    }
    suffix <- ".txt.top1000.tcrs"
    dataset <- paste(prefix, midfix, suffix, sep="_")
    return(dataset)
}

label_to_info <- function(label) {
    info_list <- label %>% strsplit(split="_") %>% unlist
    info <- paste0("Subject ", info_list[1], ", timepoint = ", info_list[3], "d")
    return(info)
}

val_df <- fread("R/data/validation.txt")
val_df[["v_gene"]] <- val_df[["v.segm"]] %>%
    sapply(discard_extra_gene_calls)

yfv_df <- val_df[val_df[["antigen.species"]] == "YellowFeverVirus", ]
cmv_df <- val_df[val_df[["antigen.species"]] == "CMV", ]

subjects <- c(
              # Timepoints in which we expect the biggest response
                 "P1_0_15",
                 "P2_0_15",
                 "Q1_0_15",
                 "Q2_0_15",
                 "S1_0_15",
                 "S2_0_15",

              # Timepoints in which we expect no response
                 "P1_0_0",
                 "P2_0_0",
                 "Q1_0_0",
                 "Q2_0_0",
                 "S1_0_0",
                 "S2_0_0",

              # Timepoints in which we expect little to no response
                 "P1_0_45",
                 "P2_0_45",
                 "Q1_0_45",
                 "Q2_0_45",
                 "S1_0_45",
                 "S2_0_45"
                )

subject_dfs <- {}
for(subject in subjects) {
    subject_df <- fread(file.path("data/yfv", label_to_dataset(subject)), header=F)
    names(subject_df) <- c("v_gene", "cdr3")
    subject_dfs[[subject]] <- subject_df
}

cluster_cols <- c("v_gene", "cdr3", "subject", "cluster")
cluster_df <- matrix(NA, nrow=0, ncol=length(cluster_cols)) %>%
    data.table %>%
    setNames(cluster_cols)

for(s in subjects) {
    dir <- file.path("output/hmm", s)
    clusters <- 1:5 %>% paste("cluster", ., sep="_")
    for(cluster in clusters) {
        tcrs <- fread(file.path(dir, cluster, "cluster_tcrs.csv"), header=T)
        v_genes <- tcrs[["tcr"]] %>% sapply(function(x) { strsplit(x, ",") %>% unlist %>% first })
        cdr3s <- tcrs[["tcr"]] %>% sapply(function(x) { strsplit(x, ",") %>% unlist %>% last })
        cluster_df <- rbind(
                            cluster_df, 
                            data.table(v_gene=v_genes, cdr3=cdr3s, subject=label_to_info(s), cluster=cluster))

    }
}

hit_cols <- c("v_gene", "cdr3", "subject", "cluster")
yfv_hits <- matrix(NA, nrow=0, ncol=length(hit_cols)) %>%
    data.table %>%
    setNames(hit_cols)
cmv_hits <- matrix(NA, nrow=0, ncol=length(hit_cols)) %>%
    data.table %>%
    setNames(hit_cols)

occ_cols <- c("v_gene", "cdr3", "subject")
yfv_occ <- matrix(NA, nrow=0, ncol=length(occ_cols)) %>%
    data.table %>%
    setNames(occ_cols)
cmv_occ <- matrix(NA, nrow=0, ncol=length(occ_cols)) %>%
    data.table %>%
    setNames(occ_cols)

for(i in 1:nrow(yfv_df)) {
    print(i)
    v_gene <- yfv_df[i, "v_gene"] %>% unlist
    cdr3 <- yfv_df[i, "cdr3"] %>% unlist
    for(j in 1:nrow(cluster_df)) {
        if(v_gene == cluster_df[j, "v_gene"] && cdr3 == cluster_df[j, "cdr3"]) {
            yfv_hits <- rbind(
                              yfv_hits,
                              cluster_df[j, ]
                             )
        }
    }

    for(s in subjects) {
        subject_df <- subject_dfs[[s]]
        for(j in 1:nrow(subject_df)) {
            if(subject_df[j, "v_gene"] == v_gene && subject_df[j, "cdr3"] == cdr3) {
                yfv_occ <- rbind(
                                 yfv_occ,
                                 data.table(
                                            v_gene=v_gene,
                                            cdr3=cdr3,
                                            subject=s
                                           )
                                )
            }
        }
    }
}

for(i in 1:nrow(cmv_df)) {
    print(i)
    v_gene <- cmv_df[i, "v_gene"] %>% unlist
    cdr3 <- cmv_df[i, "cdr3"] %>% unlist
    for(j in 1:nrow(cluster_df)) {
        if(v_gene == cluster_df[j, "v_gene"] && cdr3 == cluster_df[j, "cdr3"]) {
            cmv_hits <- rbind(
                              cmv_hits,
                              cluster_df[j, ]
                             )
        }
    }
    for(s in subjects) {
        subject_df <- subject_dfs[[s]]
        for(j in 1:nrow(subject_df)) {
            if(subject_df[j, "v_gene"] == v_gene && subject_df[j, "cdr3"] == cdr3) {
                cmv_occ <- rbind(
                                 cmv_occ,
                                 data.table(
                                            v_gene=v_gene,
                                            cdr3=cdr3,
                                            subject=s
                                           )
                                )
            }
        }
    }
}
