source("plot_utils.R")

if(!exists("fg_dat")) {
    source("cutoff.R")
}

if(!exists("full_dist_mat")) {
    full_dist_mat <- read.csv(
                              paste0(
                                     "~/sync/per_tcr/dist_matrices/",
                                     "all_subjects",
                                     ".csv"
                                    ),
                              header=FALSE
                             )
}

if(!exists("full_projection")) {
    full_projection <- get_projection(
            full_dist_mat,
            projection_method=projection_method,
            k=2
        )
}

mds_df <- build_mds_dataframe(fg_dat, radius=50.5, subjects=subjects, add_extra_metrics=TRUE)
full_dist_df <- cbind(
    full_projection,
    mds_df[!duplicated(mds_df[["tcr"]]), ]
)
names(full_dist_df)[c(1, 2)] <- c("X", "Y")

full_mds_plot <- full_dist_df %>%
    ggplot(aes(x=X, y=Y, color=relative_score, size=relative_score, shape=label)) +
    geom_point() +
    scale_colour_viridis_c() +
    scale_size(range=c(0.1, 2))
ggsave("full_mds_plot.pdf")
