library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
# source("R/plot_utils.R")

library("rjson")
CONFIG <- fromJSON(file = "config.json")
cluster_df_dir <- CONFIG["IEL_CLUSTER_OUTPUT"]

if(!exists("fg_dat")) {
    source("R/load_score_datasets.R")
}

fg_dat_sub <- fg_dat[fg_dat[["TCRDistRadius"]] == 48.5, ]
bg_dat_sub <- bg_dat[bg_dat[["TCRDistRadius"]] == 48.5, ]

df_colnames <- c("score", "tcr", "ecdf", "subject", "motif", "group")
motif_levels <- c("OT-Tremont", "OT-Revere", "OT-Ida", "N/A")

cluster_df <- matrix(NA, nrow=0, ncol=length(df_colnames)) %>%
    data.table %>%
    setNames(df_colnames)
cluster_df[["motif"]] <- factor(cluster_df[["motif"]], levels=motif_levels)

prevalence_df_colnames <- c("motif", "prevalence", "group", "subject", "count")
prevalence_df <- matrix(NA, nrow=0, ncol=length(prevalence_df_colnames)) %>%
    data.table %>%
    setNames(prevalence_df_colnames)

groups <- c("DN", "CD4", "CD8")
motifs <- c("OT-Tremont", "OT-Revere", "OT-Ida")
dn_subjects <- 1:23 %>% sapply(function(x) { paste0("DN_", x) })
cd4_subjects <- 1:23 %>% sapply(function(x) { paste0("DN_", x) })
cd8_subjects <- 1:23 %>% sapply(function(x) { paste0("DN_", x) })

for(group in groups) {
    group_subjects <- 1:23 %>% sapply(function(x) { paste(group, x, sep="_") })
    for(subject in group_subjects) {
        subject_df <- file.path(cluster_df_dir, subject, "cluster_df.csv") %>%
            data.table::fread(header=FALSE)  %>%
            setNames(c("v_gene", "cdr3", "motif"))
        subject_df[["motif"]] <- subject_df[["motif"]] %>%
            sapply(
                function(x) { 
                    if(x != "N/A") { 
                        paste("OT", x, sep="-")
                    } else {
                        "N/A"
                    }
                } 
            )
        for(motif in motifs) {
            prevalence <- mean(subject_df[["motif"]] == motif)
            counts <- sum(subject_df[["motif"]] == motif)
            prevalence_df <- prevalence_df %>% rbind(data.table(motif=motif, prevalence=prevalence, group=group, subject=subject, count=counts))
        }
        if(group == "DN") {
            subject_df <- subject_df[!duplicated(subject_df)]

            bg_loneliness_df <- bg_dat_sub[bg_dat_sub[["Subject"]] == paste0(subject, "_B"), ]
            bg_ecdf_function <- ecdf(bg_loneliness_df[["Score"]])
            bg_loneliness_df[["ecdf"]] <- bg_loneliness_df[["Score"]] %>% bg_ecdf_function

            subject_df <- cbind.data.frame(
                                score=bg_loneliness_df[["Score"]],
                                tcr=bg_loneliness_df[["TCR"]],
                                ecdf=bg_loneliness_df[["ecdf"]],
                                subject=subject,
                                motif=subject_df[["motif"]],
                                group="background"
                               )
            cluster_df <- rbind(cluster_df, subject_df)

            fg_loneliness_df <- fg_dat_sub[fg_dat_sub[["Subject"]] == paste0(subject, "_B"), ]
            fg_ecdf_function <- ecdf(fg_loneliness_df[["Score"]])
            fg_loneliness_df[["ecdf"]] <- fg_loneliness_df[["Score"]] %>% fg_ecdf_function

            subject_df <- cbind.data.frame(
                                score=fg_loneliness_df[["Score"]],
                                tcr=fg_loneliness_df[["TCR"]],
                                ecdf=fg_loneliness_df[["ecdf"]],
                                subject=subject,
                                motif=subject_df[["motif"]],
                                group="foreground"
                               )

            cluster_df <- rbind(cluster_df, subject_df)
        }

    }
}
prevalence_df[["motif"]] <- factor(prevalence_df[["motif"]], levels=motif_levels)
prevalence_df[["group"]] <- factor(prevalence_df[["group"]], levels=c("DN", "CD4", "CD8"))

outdir <- "output/cluster_iels/motif_metrics"
dir.create(outdir)

motif_colors <- brewer.pal(4, "Dark2")
names(motif_colors) <- levels(prevalence_df[["motif"]])
color_scale <- scale_colour_manual(name = "cluster", values = motif_colors)

# Rename motif to cluster for consistency in the manuscript
cluster_df[["cluster"]] <- cluster_df[["motif"]]
cluster_df[["motif"]] <- NULL
prevalence_df[["cluster"]] <- prevalence_df[["motif"]]
prevalence_df[["motif"]] <- NULL


p_ecdf <- cluster_df %>% 
    ggplot(aes(x=ecdf, y=..density.., color=cluster, linetype=cluster)) + 
    geom_freqpoly() +
    facet_wrap(vars(group)) +
    theme_minimal() +
    color_scale
ggsave(file.path(outdir, "ecdf.pdf"), width=8, height=3)

p_score <- cluster_df %>% 
    ggplot(aes(x=score, y=..density.., color=cluster, linetype=cluster)) + 
    geom_freqpoly() +
    facet_wrap(vars(group)) +
    theme_minimal() +
    color_scale
ggsave(file.path(outdir, "score.pdf"), width=8, height=3)

p_prevalence <- prevalence_df %>%
    ggplot(aes(x=prevalence, y=..density.., group=interaction(cluster, group), colour=cluster, linetype=group)) +
    geom_freqpoly(bins=10) +
    facet_wrap(vars(cluster)) +
    theme_minimal() +
    xlim(c(0, max(prevalence_df[["prevalence"]]))) +
    color_scale
ggsave(file.path(outdir, "prevalence.pdf"), width=10, height=3)

p_counts <- prevalence_df %>%
    ggplot(aes(x=count, y=..density.., group=interaction(cluster, group), colour=cluster, linetype=group)) +
    geom_freqpoly(bins=10) +
    facet_wrap(vars(cluster)) +
    theme_minimal() +
    color_scale
ggsave(file.path(outdir, "counts.pdf"), width=10, height=4)
