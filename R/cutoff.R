library(cowplot)
library(dplyr)
library(ggplot2)
library(magrittr)
library(reshape2)
library(rjson)

source("plot_utils.R")

projection_method <- "MDS"
dir.create(projection_method)

if(T) {
    dist_mats <- read.csv("dist_mat.csv", header=FALSE)
    results <- list("results"=fromJSON(file="empirical_fg_bg_nbhd_stats.json"))
    for(subject in names(results)) {
        results[[subject]][["sds"]] <- NULL
    }
    dat <- results %>% melt

    dat[, c("L3", "L1")] <- NULL
    names(dat) <- c("Score", "Group", "NeighborCount", "TCRDistRadius", "TCR", "Index", "Subject")
    dat[["NeighborCount"]] <- dat[["NeighborCount"]] %>% as.numeric
    dat[["TCRDistRadius"]] <- dat[["TCRDistRadius"]] %>% as.numeric
    dat[["Subject"]] <- dat[["Subject"]] %>% 
        sapply(gsub, pattern=".tcrs", replacement="")
    
    dat[dat$Group == "background", ]$Group <- "tmp"
    dat[dat$Group == "foreground", ]$Group <- "background"
    dat[dat$Group == "tmp", ]$Group <- "foreground"

    fg_dat <- dat[dat$Group == "foreground", ]
    bg_dat <- dat[dat$Group == "background", ]
    
    subjects <- dat$Subject %>% 
        unique
}

#rand_dict <- fromJSON(file=file.path("neighborhood_sums", "per_tcr.json"))
#rand_df <- rand_dict %>% melt
#names(rand_df) <- c("Score", "CDR3", "Gene", "TCR_ID", "TCR"F)
#rand_df[["TCR"]] <- Map(function(x, y) { paste(x, y, sep=",") }, rand_df[["Gene"]], rand_df[["CDR3"]])
#stop()

if(T) {
    dist_mats <- {}
    for(subject in subjects) {
        dist_mats[[subject]] <- read.csv(
                                         paste0(
                                                   "~/sync/per_tcr/dist_matrices/",
                                                   subject,
                                                   ".tcrs.csv"
                                                  ),
                                         header=FALSE
                                        )
    }
}

if(T) {
    mds_dats <- {}
    for(subject in subjects) {
        mds_dats[[subject]] <- get_projection(dist_mats[[subject]], 
                                              projection_method=projection_method, 
                                              k=2)
    }
}

if(FALSE) {
    xs <- mds_fit[, 1]
    ys <- mds_fit[, 2]
    scores <- fg_dat[fg_dat$TCRDistRadius == 50.5, ]$Score
    tcrs <- (dat$TCR %>% rle %$% values)
    
    left_tcr_cluster <- tcrs[xs > -70 & xs < -50 & ys > 20 & ys < 40 & scores > 1500]
    right_tcr_cluster <- tcrs[xs > 30 & xs < 60 & ys > 50 & ys < 90 & scores > 1000]
    tcr_band <- tcrs[ys > -50 & ys < -20 & scores > 1000]
}

if(TRUE) {
    tmp <- fg_dat %>% build_mds_dataframe(radius=50.5, subjects="DN_12_B")
    tcrs <- dat[dat[["Subject"]] == "DN_12_B", ][["TCR"]] %>% rle %$% values
    tcrs[tmp[tmp$x1 > -30 & tmp$x1 < 20 & tmp$x2 < -50 & tmp$score > 1800, ] %>% 
         rownames %>% 
         as.numeric]
}

fg_plots <- {}
bg_plots <- {}
radii <- fg_dat$TCRDistRadius %>% unique
for(i in 1:length(radii)) {
    fg_plots[[i]] <- plot_mds_with_scores(fg_dat, radii[i], subjects[2])
    bg_plots[[i]] <- plot_mds_with_scores(bg_dat, radii[i], subjects[2])
}
fg_mds_plots <- plot_grid(plotlist=fg_plots)
ggsave("fg_mds.pdf", width=12, height=10)
bg_mds_plots <- plot_grid(plotlist=bg_plots)
ggsave("bg_mds.pdf", width=12, height=10)

fg_mds_dat_by_subject <- fg_dat %>% 
    build_mds_dataframe(radius=50.5, subjects=subjects, add_extra_metrics=TRUE) %>%
    cbind.data.frame(group="foreground")
bg_mds_dat_by_subject <- bg_dat %>% 
    build_mds_dataframe(radius=50.5, subjects=subjects, add_extra_metrics=TRUE) %>%
    cbind.data.frame(group="background")
mds_dat_by_subject <- rbind.data.frame(
                                       fg_mds_dat_by_subject,
                                       bg_mds_dat_by_subject
                                      )

motif_metrics_dir <- "motif_metrics"
dir.create(motif_metrics_dir)

mds_dat_by_subject %>%
    ggplot(aes(x=score, y=..density.., color=label)) + 
    facet_wrap(vars(group), dir="v") +
    geom_freqpoly() +
    xlab("Average loneliness") +
    theme_minimal()
ggsave(file.path(motif_metrics_dir, "loneliness_by_motif.pdf"))

mds_dat_by_subject %>%
    ggplot(aes(x=relative_score, y=..density.., color=label)) + 
    facet_wrap(vars(group), dir="v") +
    geom_freqpoly() +
    xlab("Average relative loneliness") +
    theme_minimal()
ggsave(file.path(motif_metrics_dir, "relative_loneliness_by_motif.pdf"))

mds_dat_by_subject %>%
    ggplot(aes(x=ecdf, y=..density.., color=label)) + 
    facet_wrap(vars(group)) +
    geom_freqpoly() +
    xlab("ECDF") +
    theme_minimal()
ggsave(file.path(motif_metrics_dir, "ecdf_by_motif.pdf"))

prevalence_dat <- mds_dat_by_subject[, c("label", "prevalence")]
prevalence_dat[!duplicated(prevalence_dat), ] %>%
    ggplot(aes(x=prevalence, y=..density.., color=label)) + 
    geom_freqpoly() +
    xlab("Prevalence") +
    theme_minimal()
ggsave(file.path(motif_metrics_dir, "motif_prevalence_fg.pdf"))

snapshot_dir <- file.path("snapshots", projection_method)
dir.create(snapshot_dir)
for(subject in subjects) {
    fg_snapshot_dat  <- fg_dat %>%
       build_mds_dataframe(radius=50.5, subjects=subject) %>%
       cbind.data.frame(Group="foreground")
    bg_snapshot_dat  <- bg_dat %>%
       build_mds_dataframe(radius=50.5, subject=subject) %>%
       cbind.data.frame(Group="background")
    snapshot_dat <- rbind.data.frame(fg_snapshot_dat, bg_snapshot_dat)
    snapshot_dat[["label"]] <- factor(snapshot_dat[["label"]], levels=c("N/A", "Revere", "Tremont", "Ida", "X"))
    snapshot_dat %>%
        ggplot(aes(x=x1, y=x2, color=score)) + geom_point(aes(shape=label)) +
           scale_colour_viridis_c()  +
           theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
                 panel.background = element_blank()) +
           facet_wrap(vars(Group))
    ggsave(file.path(snapshot_dir, paste0(subject, ".pdf")), width=10, height=6)
}

fg_mds_dat_by_subject <- fg_dat %>% build_mds_dataframe(radius=50.5, subjects=subjects)
bg_mds_dat_by_subject <- bg_dat %>% build_mds_dataframe(radius=50.5, subjects=subjects)
by_subject_scores <- as.numeric(c(fg_mds_dat_by_subject$score, bg_mds_dat_by_subject$score))

bg_dat %>%
    build_mds_dataframe(radius=50.5, subjects=subjects) %>%
    ggplot(aes(x=x1, y=x2, color=score)) + 
       geom_point(size=0.5) +
       scale_colour_viridis_c(breaks=by_subject_scores) +
       theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
             panel.background = element_blank()) +
       facet_wrap(vars(subject), scales="free")
ggsave("bg_snapshot_by_subject.pdf", width=12, height=10)

plot_dir <- "~/sync/per_tcr/radius_analysis"

p <- dat %>% ggplot(aes(x=Score, y=TCRDistRadius, color=NeighborCount)) + 
    geom_point(size=0.5) + 
    facet_wrap(vars(Group)) + 
    xlab("Mean per-tcr loneliness") + 
    ylab("TCRdist radius")
ggsave(file.path(plot_dir, "score_vs_radius.pdf"), width=8, height=6)

p2 <- dat %>% ggplot(aes(x=Score, y=NeighborCount, colour=Group)) + 
    geom_point(size=0.5) +
    xlab("Mean per-tcr loneliness") +
    ylab("Neighbor count")
ggsave(file.path(plot_dir, "score_vs_neighbors.pdf"), width=8, height=6)

p2_tcr <- dat[dat$TCR %in% 
              (dat$TCR %>% unique %>% sample(20)), ] %>% 
    ggplot(aes(x=Score, y=NeighborCount, colour=Group)) + 
    geom_point(size=0.5) +
    xlab("Mean per-tcr loneliness") +
    ylab("Neighbor count") +
    facet_wrap(vars(TCR)) +
    theme(legend.position="none")
ggsave(file.path(plot_dir, "score_vs_neighbors_by_tcr.pdf"), width=12, height=8)

p2_facet <- dat %>% ggplot(aes(x=Score, y=NeighborCount, colour=Group)) + 
    geom_point(size=0.5) + 
    facet_wrap(vars(TCRDistRadius)) +
    xlab("Mean per-tcr loneliness") +
    ylab("Neighbr count")
ggsave(file.path(plot_dir, "score_vs_neighbors_by_cutoff.pdf"), width=12, height=8)

p2_facet_scaled <- dat %>% ggplot(aes(x=Score, y=NeighborCount, colour=Group)) + 
    geom_point(size=0.5) + 
    facet_wrap(vars(TCRDistRadius), scales="free") +
    xlab("Mean per-tcr loneliness") +
    ylab("Neighbr count")
ggsave(file.path(plot_dir, "score_vs_neighbors_by_cutoff_scaled.pdf"), width=12, height=8)

p3 <- dat %>% ggplot(aes(x=TCRDistRadius, y=NeighborCount)) + 
    geom_point(size=0.1) +
    xlab("TCRdist radius") +
    ylab("Neighbor count")
ggsave(file.path(plot_dir, "neighbors_vs_radius.pdf"), width=8, height=6)

p4 <- dat %>% ggplot(aes(x=Score, colour=Group)) + 
    stat_ecdf() + 
    facet_wrap(vars(TCRDistRadius)) +
    xlab("Mean per-tcr loneliness") +
    ylab("ECDF")
ggsave(file.path(plot_dir, "ecdfs_by_radius.pdf"))



