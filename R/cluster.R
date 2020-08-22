library(data.table)
library(dplyr)
library(ggplot2)

cluster_df_dir <- "output/iel_clusters"

if(!exists("fg_dat")) {
    source("R/load_score_datasets.R")
}

fg_dat <- fg_dat[fg_dat[["Group"]] == "foreground" & fg_dat[["TCRDistRadius"]] == 50.5, ]

df_colnames <- c("score", "tcr", "ecdf", "subject", "motif")
cluster_df <- matrix(NA, nrow=0, ncol=length(df_colnames)) %>%
    data.table %>%
    setNames(df_colnames)

prevalence_df <- matrix(NA, nrow=0, ncol=3) %>%
    data.table %>%
    setNames(c("motif", "prevalence", "group"))

groups <- c("DN", "CD4", "CD8")
motifs <- c("Ida", "Revere", "Tremont")
dn_subjects <- 1:23 %>% sapply(function(x) { paste0("DN_", x) })
cd4_subjects <- 1:23 %>% sapply(function(x) { paste0("DN_", x) })
cd8_subjects <- 1:23 %>% sapply(function(x) { paste0("DN_", x) })

for(group in groups) {
    group_subjects <- 1:23 %>% sapply(function(x) { paste(group, x, sep="_") })
    for(subject in group_subjects) {
        subject_df <- file.path(cluster_df_dir, subject, "cluster_df.csv") %>%
            data.table::fread(header=FALSE)  %>%
            setNames(c("v_gene", "cdr3", "motif"))
        if(group == "DN") {
            subject_df <- subject_df[!duplicated(subject_df)]

            loneliness_df <- fg_dat[fg_dat[["Subject"]] == paste0(subject, "_B"), ]
            ecdf_function <- ecdf(loneliness_df[["Score"]])
            loneliness_df[["ecdf"]] <- loneliness_df[["Score"]] %>% ecdf_function

            subject_df <- cbind.data.frame(
                                score=loneliness_df[["Score"]],
                                tcr=loneliness_df[["TCR"]],
                                ecdf=loneliness_df[["ecdf"]],
                                subject=subject,
                                motif=subject_df[["motif"]]
                               )

            cluster_df <- rbind(cluster_df, subject_df)
        }

        for(motif in motifs) {
            prevalence <- mean(subject_df[["motif"]] == motif)
            prevalence_df <- prevalence_df %>%
                rbind(data.table(motif=motif, prevalence=prevalence, group=group))
        }
    }
}

p_ecdf <- cluster_df %>% 
    ggplot(aes(x=ecdf, y=..density.., color=motif)) + 
    geom_freqpoly() +
    theme_minimal()

p_score <- cluster_df %>% 
    ggplot(aes(x=score, y=..density.., color=motif)) + 
    geom_freqpoly() +
    theme_minimal()

p_prevalence <- prevalence_df %>%
    ggplot(aes(x=prevalence, y=..density.., group=interaction(motif, group), colour=motif, linetype=group)) +
    geom_freqpoly(bins=15) +
    facet_wrap(vars(motif)) +
    theme_minimal()
