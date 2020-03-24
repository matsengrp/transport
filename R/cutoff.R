library(cowplot)
library(dplyr)
library(ggplot2)
library(magrittr)
library(reshape2)
library(rjson)

dist_mat <- read.csv("dist_mat.csv", header=FALSE)
results <- fromJSON(file="empirical_fg_bg_nbhd_stats.json")
results[["sds"]] <- NULL
dat <- results %>% melt
dat[["L1"]] <- NULL
names(dat) <- c("Score", "Group", "NeighborCount", "TCRDistRadius", "TCR", "Index")
dat[["NeighborCount"]] <- dat[["NeighborCount"]] %>% as.numeric
dat[["TCRDistRadius"]] <- dat[["TCRDistRadius"]] %>% as.numeric

dat[dat$Group == "background", ]$Group <- "tmp"
dat[dat$Group == "foreground", ]$Group <- "background"
dat[dat$Group == "tmp", ]$Group <- "foreground"

# Switch nomenclature from Phil's analysis
fg_dat <- dat[dat$Group == "foreground", ]
bg_dat <- dat[dat$Group == "background", ]

build_mds_dataframe <- function(ref_dat, radius) {
    mds_scores <- ref_dat[ref_dat$TCRDistRadius == radius, ]$Score
    mds_dat <- data.frame(x1=mds_fit[, 1], x2=mds_fit[, 2], score=mds_scores)
    return(mds_dat)
}

plot_mds_with_scores <- function(ref_dat, radius) {
    mds_dat <- build_mds_dataframe(ref_dat, radius)
    p <- ggplot(mds_dat, aes(x=x1, y=x2, color=score)) + 
        geom_point(size=0.5) +
        scale_colour_viridis_c(breaks = as.numeric(mds_dat$score)) +
        theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(),
              panel.background = element_blank()) +
        ggtitle(paste("Radius = ", radius))
    return(p)
}

mds_fit <- cmdscale(dist_mat, k=2)
xs <- mds_fit[, 1]
ys <- mds_fit[, 2]
scores <- fg_dat[fg_dat$TCRDistRadius == 50.5, ]$Score
tcrs <- (dat$TCR %>% rle %$% values)

left_tcr_cluster <- tcrs[xs > -70 & xs < -50 & ys > 20 & ys < 40 & scores > 1500]
right_tcr_cluster <- tcrs[xs > 30 & xs < 60 & ys > 50 & ys < 90 & scores > 1000]
tcr_band <- tcrs[ys > -50 & ys < -20 & scores > 1000]

fg_plots <- {}
bg_plots <- {}
radii <- fg_dat$TCRDistRadius %>% unique
for(i in 1:length(radii)) {
    fg_plots[[i]] <- plot_mds_with_scores(fg_dat, radii[i])
    bg_plots[[i]] <- plot_mds_with_scores(bg_dat, radii[i])
}
fg_mds_plots <- plot_grid(plotlist=fg_plots)
ggsave("fg_mds.pdf", width=12, height=10)
bg_mds_plots <- plot_grid(plotlist=bg_plots)
ggsave("bg_mds.pdf", width=12, height=10)

fg_snapshot_dat  <- fg_dat %>%
   build_mds_dataframe(radius=50.5) %>%
   cbind.data.frame(Group="foreground")
bg_snapshot_dat  <- bg_dat %>%
   build_mds_dataframe(radius=50.5) %>%
   cbind.data.frame(Group="background")
snapshot_dat <- rbind.data.frame(fg_snapshot_dat, bg_snapshot_dat)
snapshot_dat %>%
    ggplot(aes(x=x1, y=x2, color=score)) + geom_point() +
       scale_colour_viridis_c()  +
       theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
             panel.background = element_blank()) +
       facet_wrap(vars(Group))
ggsave("snapshot.pdf", width=10, height=6)

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


