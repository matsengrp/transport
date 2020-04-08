library(cowplot)
library(dplyr)
library(ggplot2)
library(magrittr)
library(reshape2)
library(rjson)

if(FALSE) {
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
    dist_mats <- {}
    mds_dats <- {}
    for(subject in subjects) {
        dist_mats[[subject]] <- read.csv(
                                         paste0(
                                                   "~/sync/per_tcr/dist_matrices/",
                                                   subject,
                                                   ".tcrs.csv"
                                                  ),
                                         header=FALSE
                                        )
        mds_dats[[subject]] <- cmdscale(dist_mats[[subject]], k=2)
    }
}


build_mds_dataframe <- function(ref_dat, radius, subjects) {
    full_dat <- matrix(NA, nrow=0, ncol=4) %>% 
        data.frame %>%
        setNames(c("x1", "x2", "subject", "score"))
    for(subject in subjects) {
        subject_radius_dat <- ref_dat[ref_dat[["TCRDistRadius"]] == radius & 
                                      ref_dat[["Subject"]] == subject, ]
        subject_dat <- data.frame(
                                  x1=mds_dats[[subject]][, 1], 
                                  x2=mds_dats[[subject]][, 2],
                                  subject=subject,
                                  score=subject_radius_dat[["Score"]],
                                  label=subject_radius_dat[["TCR"]] %>% sapply(get_motif_label)
                                 )
        full_dat <- rbind.data.frame(full_dat, subject_dat) 
    }
    return(full_dat)
}


plot_mds_with_scores <- function(ref_dat, radius, subject) {
    mds_dat <- build_mds_dataframe(ref_dat, radius, subject)
    p <- ggplot(mds_dat, aes(x=x1, y=x2, color=score)) + 
        geom_point(size=0.5) +
        scale_colour_viridis_c(breaks = as.numeric(mds_dat$score)) +
        theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(),
              panel.background = element_blank()) +
        ggtitle(paste("Radius = ", radius))
    return(p)
}

extract_tcr_info <- function(tcr) {
    return(tcr %>% strsplit(split=",") %>% first)
}

has_revere_motif <- function(tcr) {
    tcr_info <- tcr %>% 
        extract_tcr_info
    gene <- tcr_info[1]
    cdr3 <- tcr_info[2]
    revere_cdr3 <- "GT[VI]SNERLFF"
    return(grepl(revere_cdr3, cdr3))
}

has_tremont_motif <- function(tcr) {
    tcr_info <- tcr %>%
        extract_tcr_info
    gene <- tcr_info[1]
    cdr3 <- tcr_info[2]
    tremont_gene <- "TRBV16"
    tremont_cdr3 <- "DWG"
    return(grepl(tremont_gene, gene) && grepl(tremont_cdr3, cdr3))
}

has_ida_motif <- function(tcr) {
    tcr_info <- tcr %>%
        extract_tcr_info
    gene <- tcr_info[1]
    cdr3 <- tcr_info[2]
    ida_gene <- "TRBV(16|12)"
    ida_cdr3 <- "Y*A*EQ[YF]F"
    return(grepl(ida_gene, gene) && grepl(ida_cdr3, cdr3))
}

has_motif_x <- function(tcr) {
    tcr_info <- tcr %>%
        extract_tcr_info
    gene <- tcr_info[1]
    cdr3 <- tcr_info[2]
    x_gene <- "TRBV(12|13)"
    x_cdr3  <- "[AE]E*[TR]*[LQ][FY]F"
    return(grepl(x_gene, gene) && grepl(x_cdr3, cdr3))
}

get_motif_label <- function(tcr) {
    if(tcr %>% has_revere_motif) {
        label <- "Revere"
    } else if(tcr %>% has_tremont_motif) {
        label <- "Tremont"
    } else if(tcr %>% has_ida_motif) {
        label <- "Ida"
    } else if(tcr %>% has_motif_x) {
        label <- "X"
    } else {
        label <- "N/A"
    }
    return(label)
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

extract_tcrs_from_mds_cluster <- function(
                                          ref_dat,
                                          subject,
                                          xrange,
                                          yrange,
                                          score_threshold=1800,
                                          radius=50.5
                                         ) {
   tmp <- ref_dat %>%
       build_mds_dataframe(radius=radius, subjects=subject)
   filtered_dat <- tmp[tmp[["x1"]] > xrange[1] &
                       tmp[["x1"]] < xrange[2] &
                       tmp[["x2"]] > yrange[1] &
                       tmp[["x2"]] < yrange[2] &
                       tmp[["score"]] > score_threshold,
                      ]
   tcrs <- dat[dat[["Subject"]] == subject, ][["TCR"]] %>% 
       rle %$%
       values
   cluster_tcrs <- tcrs[filtered_dat %>% rownames %>% as.numeric]
   return(cluster_tcrs)
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
    ggsave(paste0("snapshots/", subject, ".pdf"), width=10, height=6)
}

fg_mds_dat_by_subject <- fg_dat %>% build_mds_dataframe(radius=50.5, subjects=subjects)
bg_mds_dat_by_subject <- bg_dat %>% build_mds_dataframe(radius=50.5, subjects=subjects)
by_subject_scores <- as.numeric(c(fg_mds_dat_by_subject$score, bg_mds_dat_by_subject$score))

fg_mds_dat_by_subject %>%
    ggplot(aes(x=x1, y=x2, color=score)) + 
       geom_point(size=0.5) +
       scale_colour_viridis_c(breaks=by_subject_scores) +
       theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
             panel.background = element_blank()) +
       facet_wrap(vars(subject), scales="free")
ggsave("fg_snapshot_by_subject.pdf", width=12, height=10)
    
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



