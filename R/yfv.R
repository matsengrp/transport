library(cowplot)
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
yfv_df[["v_gene_no_allele"]] <- yfv_df[["v_gene"]] %>%
    sapply(gsub, pattern="\\*\\d*", replacement="")
yfv_df[["tcr_no_allele"]] <- paste(yfv_df[["v_gene_no_allele"]],
                                       yfv_df[["cdr3"]],
                                       sep=","
                                      )
cmv_df <- val_df[val_df[["antigen.species"]] == "CMV", ]
cmv_df[["v_gene_no_allele"]] <- cmv_df[["v_gene"]] %>%
    sapply(gsub, pattern="\\*\\d*", replacement="")
cmv_df[["tcr_no_allele"]] <- paste(cmv_df[["v_gene_no_allele"]],
                                       cmv_df[["cdr3"]],
                                       sep=","
                                      )

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
    subject_df[["v_gene_no_allele"]] <- subject_df[["v_gene"]] %>%
        sapply(gsub, pattern="\\*\\d*", replacement="")
    subject_df[["tcr_no_allele"]] <- paste(subject_df[["v_gene_no_allele"]],
                                           subject_df[["cdr3"]],
                                           sep=","
                                          )
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
cluster_df[["v_gene_no_allele"]] <- cluster_df[["v_gene"]] %>%
    sapply(gsub, pattern="\\*\\d*", replacement="")
cluster_df[["tcr"]] <- paste(cluster_df[["v_gene_no_allele"]], cluster_df[["cdr3"]], sep=",")

common_hits_colnames <- c("subject", "timepoint", "cluster", "hit_rate")
common_hits <- matrix(NA, nrow=0, ncol=length(common_hits_colnames)) %>%
    data.table %>%
    setNames(common_hits_colnames)
common_hits_without_subject <- common_hits
common_hits_without_subject[["subject"]] <- NULL
common_hits_without_cluster <- common_hits
common_hits_without_cluster[["cluster"]] <- NULL
common_hits_without_cluster$"cluster_cutoff" <- NA
common_hits_without_subject$"agg_hit_rate" <- NA
cluster_df[["cluster_int"]] <- cluster_df[["cluster"]] %>%
    sapply(function(x) { strsplit(x, split="_") %>% unlist %>% last %>% as.numeric })
for(tmp_timepoint in cluster_df[["timepoint"]] %>% unique) {
    for(tmp_cluster in clusters) {
        sub_df <- cluster_df[
                             cluster_df[["timepoint"]] == tmp_timepoint &
                             cluster_df[["cluster"]] == tmp_cluster, ]
        agg_sub_df <- cluster_df[
                             cluster_df[["timepoint"]] == tmp_timepoint &
                             cluster_df[["cluster_int"]] <= strsplit(tmp_cluster, split="_") %>% unlist %>% last %>% as.numeric, ]
        hits <- {}
        agg_hits <- {}
        for(tmp_subject in cluster_df[["subject"]] %>% unique) {
            hits <- c(hits, sub_df[sub_df[["subject"]] == tmp_subject, ][["tcr"]] %in% pog_df[pog_df[["donor"]] == tmp_subject, ][["tcr"]])
            agg_hits <- c(agg_hits, agg_sub_df[agg_sub_df[["subject"]] == tmp_subject, ][["tcr"]] %in% pog_df[pog_df[["donor"]] == tmp_subject, ][["tcr"]])
        }
        hit_rate <- hits %>% mean
        agg_hit_rate <- agg_hits %>% mean
        common_hits_without_subject <- rbind(
                             common_hits_without_subject,
                             data.table(
                                        timepoint=paste0("0d vs ", tmp_timepoint, "d"),
                                        cluster=tmp_cluster,
                                        hit_rate=hit_rate,
                                        agg_hit_rate=agg_hit_rate
                                       )
                             )
    }
    for(tmp_subject in cluster_df[["subject"]] %>% unique) {
        for(cluster_cutoff in c(2, 10)) {
            sub_df <- cluster_df[
                                 cluster_df[["cluster_int"]] <= cluster_cutoff &
                                 cluster_df[["timepoint"]] == tmp_timepoint &
                                 cluster_df[["subject"]] == tmp_subject, ]
            hit_rate <- sub_df[["tcr"]] %in% pog_df[pog_df[["donor"]] == tmp_subject, ][["tcr"]] %>% mean
            common_hits_without_cluster <- rbind(
                                 common_hits_without_cluster,
                                 data.table(
                                            subject=tmp_subject,
                                            timepoint=paste0("0d vs ", tmp_timepoint, "d"),
                                            hit_rate=hit_rate,
                                            cluster_cutoff=cluster_cutoff
                                           )
                                 )
        }
    }
}
common_hits[["cluster"]] <- common_hits[["cluster"]] %>%
    sapply(gsub, pattern="cluster_", replacement="") %>%
    factor(levels=1:length(clusters))
common_hits_without_subject[["cluster"]] <- common_hits_without_subject[["cluster"]] %>%
    sapply(gsub, pattern="cluster_", replacement="") %>%
    factor(levels=1:length(clusters))

plot_dir <- "output/yfv"
p_hits_no_cluster_2 <- common_hits_without_cluster[cluster_cutoff == 2, ] %>%
    ggplot(aes(x=subject, y=hit_rate, fill=subject)) +
    geom_bar(stat="identity") +
    ylab("subject") +
    ylim(0, 1) +
    facet_wrap(vars(timepoint)) +
    theme_minimal_hgrid()
ggsave(file.path(plot_dir, "hits_by_subject_2.pdf"), width=7, height=2)
p_hits_no_cluster_10 <- common_hits_without_cluster[cluster_cutoff == 10, ] %>%
    ggplot(aes(x=subject, y=hit_rate, fill=subject)) +
    geom_bar(stat="identity") +
    ylab("subject") +
    ylim(0, 1) +
    facet_wrap(vars(timepoint)) +
    theme_minimal_hgrid()
ggsave(file.path(plot_dir, "hits_by_subject_10.pdf"), width=7, height=2)
p_hits_no_subject <- common_hits_without_subject %>%
    ggplot(aes(x=cluster, y=hit_rate)) +
    geom_bar(stat="identity") +
    xlab("cluster") +
    ylab("hit rate") +
    ylim(0, 1) +
    facet_wrap(vars(timepoint)) +
    theme_minimal_hgrid()
ggsave(file.path(plot_dir, "hits_by_cluster.pdf"), width=7, height=2)
p_agg_hits_no_subject <- common_hits_without_subject %>%
    ggplot(aes(x=cluster, y=agg_hit_rate)) +
    geom_bar(stat="identity") +
    #geom_point() +
    #geom_line(aes(x=as.numeric(cluster), y=agg_hit_rate)) +
    xlab("cluster") +
    ylab("aggregate hit rate") +
    ylim(0, 1) +
    facet_wrap(vars(timepoint)) +
    theme_minimal_hgrid()
ggsave(file.path(plot_dir, "agg_hits_by_cluster.pdf"), width=7, height=2)
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

hit_filename <- "output/hmm/yfv_hits.csv"
hit_df <- fread(hit_filename)
hit_df[['timepoint']] <- hit_df[['comparison']] %>%
    sapply(label_to_timepoint)
hit_df[['subject']] <- hit_df[['comparison']] %>%
    sapply(label_to_subject)
hit_df[['cluster']] <- hit_df[['cluster']] %>% 
    sapply(function(x) { x %>% strsplit(split="_") %>% unlist %>% last}) %>%
    as.numeric
hit_df_top_3 <- hit_df[hit_df[['timepoint']] == "15" & hit_df[["cluster"]] <= 3, ]
hit_df_top_10 <- hit_df[hit_df[['timepoint']] == "15" & hit_df[["cluster"]] <= 10, ]

occ_cols <- c("v_gene", "cdr3", "tcr", "subject")
yfv_occ <- matrix(NA, nrow=0, ncol=length(occ_cols)) %>%
    data.table %>%
    setNames(occ_cols)
cmv_occ <- matrix(NA, nrow=0, ncol=length(occ_cols)) %>%
    data.table %>%
    setNames(occ_cols)

cluster_df[["tcr"]] <- paste(cluster_df[["v_gene"]], cluster_df[["cdr3"]], sep=",")

yfv_hits <- cluster_df[tcr %in% yfv_df[['tcr']] & timepoint == "15", ]
cmv_hits <- cluster_df[tcr %in% cmv_df[['tcr']] & timepoint == "15", ]
