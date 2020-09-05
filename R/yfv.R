library(data.table)
library(dplyr)
library(ggplot2)

info_to_subject <- function(info) {
    split_info <- info %>% strsplit(split=" ") %>% unlist
    return(split_info[2] %>% gsub(pattern=",", replacement=""))
}

discard_extra_gene_calls <- function(genes) {
    top_gene <- genes %>% 
        sapply(strsplit, split=",") %>% 
        unlist %>% 
        first
    return(top_gene)
}

label_to_subject <- function(label) {
    info_list <- label %>% strsplit(split="_") %>% unlist
    return(info_list %>% first)
}

label_to_timepoint <- function(label) {
    info_list <- label %>% strsplit(split="_") %>% unlist
    return(info_list %>% last)
}

label_to_dataset <- function(label) {
    info_list <- label %>% strsplit(split="_") %>% unlist
    prefix <- info_list[1]
    timepoint <- info_list[3]
    if(timepoint == 0) {
        midfix <- "0_F2"
    } else if(timepoint == -7) {
        midfix <- "pre0_F1"
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
val_df[["tcr"]] <- paste(val_df[["v_gene"]], val_df[["cdr3"]], sep=",")

yfv_df <- val_df[val_df[["antigen.species"]] == "YellowFeverVirus", ]
cmv_df <- val_df[val_df[["antigen.species"]] == "CMV", ]

pog_df <- fread("R/data/pogorelyy_tcrs.txt")
pog_df[["tcr"]] <- paste(pog_df[["bestVGene"]], pog_df[["CDR3.amino.acid.sequence"]], sep=",")


subjects <- c(
              # Timepoints in which we expect the biggest response
                 "P1_0_15",
                 "P2_0_15",
                 "Q1_0_15",
                 "Q2_0_15",
                 "S1_0_15",
                 "S2_0_15",

              # Timepoints in which we expect no response
                 #"P1_0_0",
                 #"P2_0_0",
                 #"Q1_0_0",
                 #"Q2_0_0",
                 #"S1_0_0",
                 #"S2_0_0",

                 "P1_0_-7",
                 "P2_0_-7",
                 "Q1_0_-7",
                 "Q2_0_-7",
                 "S1_0_-7",
                 "S2_0_-7",

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
    subject_df[["tcr"]] <- paste(subject_df[["v_gene"]], subject_df[["cdr3"]], sep=",")
    subject_dfs[[subject]] <- subject_df
}

cluster_cols <- c("v_gene", "cdr3", "subject", "timepoint", "cluster")
cluster_df <- matrix(NA, nrow=0, ncol=length(cluster_cols)) %>%
    data.table %>%
    setNames(cluster_cols)

clusters <- 1:10 %>% paste("cluster", ., sep="_")
for(s in subjects) {
    dir <- file.path("output/hmm", s)
    for(cluster in clusters) {
        tcrs <- fread(file.path(dir, cluster, "cluster_tcrs.csv"), header=T)
        if(nrow(tcrs) > 3) {
        v_genes <- tcrs[["tcr"]] %>% sapply(function(x) { strsplit(x, ",") %>% unlist %>% first })
        cdr3s <- tcrs[["tcr"]] %>% sapply(function(x) { strsplit(x, ",") %>% unlist %>% last })
        cluster_df <- rbind(
                            cluster_df, 
                            data.table(
                                       v_gene=v_genes,
                                       cdr3=cdr3s,
                                       subject=label_to_subject(s),
                                       timepoint=label_to_timepoint(s),
                                       cluster=cluster
                                      )
                            )
        }

    }
}
cluster_df[["v_gene_no_allele"]] <- cluster_df[["v_gene"]] %>%
    sapply(gsub, pattern="\\*\\d*", replacement="")
cluster_df[["tcr"]] <- paste(cluster_df[["v_gene_no_allele"]], cluster_df[["cdr3"]], sep=",")

common_hits_colnames <- c("subject", "timepoint", "cluster", "hit_rate")
common_hits <- matrix(NA, nrow=0, ncol=length(common_hits_colnames)) %>%
    data.table %>%
    setNames(common_hits_colnames)
for(tmp_subject in cluster_df[["subject"]] %>% unique) {
    for(tmp_timepoint in cluster_df[["timepoint"]] %>% unique) {
        for(tmp_cluster in clusters) {
            sub_df <- cluster_df[cluster_df[["subject"]] == tmp_subject & 
                                 cluster_df[["timepoint"]] == tmp_timepoint &
                                 cluster_df[["cluster"]] == tmp_cluster, ]
            hit_rate <- sub_df[["tcr"]] %in% pog_df[pog_df[["donor"]] == tmp_subject, ][["tcr"]] %>% mean
            common_hits <- rbind(
                                 common_hits,
                                 data.table(subject=tmp_subject,
                                            timepoint=tmp_timepoint,
                                            cluster=tmp_cluster,
                                            hit_rate=hit_rate
                                           )
                                 )
        }
    }
}
common_hits[["cluster"]] <- common_hits[["cluster"]] %>%
    sapply(gsub, pattern="cluster_", replacement="") %>%
    factor(levels=1:length(clusters))
p_hits <- common_hits %>%
    ggplot(aes(x=cluster, y=hit_rate)) +
    geom_boxplot(outlier.size=-1) +
    facet_wrap(vars(timepoint)) + 
    geom_point(aes(color=subject, shape=subject))

hit_cols <- c("tcr", "subject", "cluster")
yfv_hits <- matrix(NA, nrow=0, ncol=length(hit_cols)) %>%
    data.table %>%
    setNames(hit_cols)
cmv_hits <- matrix(NA, nrow=0, ncol=length(hit_cols)) %>%
    data.table %>%
    setNames(hit_cols)

occ_cols <- c("v_gene", "cdr3", "tcr", "subject")
yfv_occ <- matrix(NA, nrow=0, ncol=length(occ_cols)) %>%
    data.table %>%
    setNames(occ_cols)
cmv_occ <- matrix(NA, nrow=0, ncol=length(occ_cols)) %>%
    data.table %>%
    setNames(occ_cols)

yfv_hit_locs <- which(yfv_df[["tcr"]] %in% cluster_df[["tcr"]])
yfv_hits <- yfv_df[which(yfv_df[["tcr"]] %in% cluster_df[["tcr"]]), "tcr"]
yfv_hits[["subject"]] <- yfv_hits[["tcr"]] %>%
    sapply(function(x) { cluster_df[tcr == x, "subject"] }) %>%
    unlist
yfv_hits[["timepoint"]] <- yfv_hits[["tcr"]] %>%
    sapply(function(x) { cluster_df[tcr == x, "timepoint"] }) %>%
    unlist
yfv_hits[["cluster"]] <- yfv_hits[["tcr"]] %>%
    sapply(function(x) { cluster_df[tcr == x, "timepoint"] }) %>%
    unlist


stop()
#for(i in 1:nrow(yfv_df)) {
#    print(i)
#    v_gene <- yfv_df[i, "v_gene"] %>% unlist
#    cdr3 <- yfv_df[i, "cdr3"] %>% unlist
#    for(j in 1:nrow(cluster_df)) {
#        if(v_gene == cluster_df[j, "v_gene"] && cdr3 == cluster_df[j, "cdr3"]) {
#            yfv_hits <- rbind(
#                              yfv_hits,
#                              cluster_df[j, ]
#                             )
#        }
#    }
#
    if(T) {
        for(s in subjects) {
            subject_df <- subject_dfs[[s]]
            yfv_occ <- rbind(
                             yfv_occ,
                             cbind(
                                   yfv_df[which(yfv_df[["tcr"]] %in% subject_df[["tcr"]]), ],
                                   subject=label_to_subject(s)
                                  )
                            )
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
    if(FALSE) {
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
}
