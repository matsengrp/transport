library(data.table)
library(dplyr)
library(ggplot2)

if(!exists("full_dist_mat")) {
    dist_mat_dir <- "output/dist_matrices"
    full_dist_mat <- read.csv(
        file.path(
            dist_mat_dir,
            paste0(
                "full_dn",
                ".csv"
            )
        ),
        header=FALSE
    )
}

if(!exists("full_projection")) {
    full_projection <- get_projection(
        full_dist_mat,
        projection_method="MDS",
        k=2
    )
}

full_projection <- full_projection %>%
    as.data.frame %>%
    setNames(c("X", "Y"))
combined_df <- fread("output/csv/three_cluster.csv")

full_mds_plot <- cbind.data.frame(full_projection, combined_df) %>%
    ggplot(aes(x=X, y=Y, color=score, size=ifelse(cluster, 2, 1), shape=as.factor(cluster))) +
    geom_point() +
    scale_colour_viridis_c() +
    scale_size(range=c(1, 2))
snapshot_dir = "output/snapshots"
snapshot_dir %>% dir.create
ggsave(file.path(snapshot_dir, "all_dn_mds_plot.pdf"), width=20, height=20)
