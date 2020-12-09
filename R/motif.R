library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(rjson)
library(segmented)
library(viridis)

source("R/plot_utils.R")

get_breakpoint_from_model <- function(seg_fit) {
    return(seg_fit[["psi"]][2])
}

csv_dir <- "output/csv"
json_dir <- "output/json"

motif_dir <- "output/motif"
e_value_motif_dir <- "output/motif/e_values"
dir.create(motif_dir)
dir.create(e_value_motif_dir)


motif_dat <- fread(file.path(csv_dir, "motif.csv"))
e_value_dat <- fread(file.path(csv_dir, "ida_e_value.csv"))

score_results <- list("results"=fromJSON(file=file.path(json_dir, "empirical_fg_bg_nbhd_stats.json"))) 

breakpoints <- {}
pdf(file.path(motif_dir, "enrichment_by_radius.pdf"), width=10, height=12)
par(mfrow=c(6, 4))
motif_dat[["subject"]] <- motif_dat[["subject"]] %>% sapply(gsub, pattern=".tcrs", replacement="")
# Enforce an ordering for the subjects for plotting
subjects <- factor(subjects, levels=1:23 %>% paste("DN", ., "B", sep="_"))
for(tmp_subject in subjects[order(subjects)]) {
    d_sub <- motif_dat[motif_dat$subject == tmp_subject, ]
    sample_size <- length(score_results[['results']][[paste0(tmp_subject, '.tcrs')]][['48.5']])
    print(sample_size)
    if(sample_size) {
        print(tmp_subject)
        lm_fit <- lm(annulus_enrichment ~ radius, data=d_sub)
        seg_fit <- segmented(lm_fit, psi=list(radius=40.5))
        breakpoint <- seg_fit %>% get_breakpoint_from_model
        breakpoints <- c(breakpoints, breakpoint)
        cutoff_radius <- radii[which(radii <= breakpoint) %>% max]
        cluster_tcrs[[tmp_subject]] <- json_object[[tmp_subject]][[toString(cutoff_radius)]][["tcrs"]]
        xs <- seq(0, max(d_sub$radius), length.out=1000)
        plot(d_sub$annulus_enrichment ~ d_sub$radius, pch=19, xlab="Radius", ylab="Mean loneliness", main=gsub(gsub(tmp_subject, pattern="DN_", replacement="Subject "), pattern="_B", replacement=""))
        lines(predict(seg_fit, newdata=data.frame(radius=xs)) ~ xs, col="red", lty=ifelse(sample_size > 200, 1, 2))
    }
}
dev.off()

stop()

pdf(file.path(motif_dir, "enrichment_by_cluster_size.pdf"), width=12, height=8)
par(mfrow=c(4, 4))
for(tmp_subject in subjects) {
    d_sub <- motif_dat[motif_dat$subject == tmp_subject, ]
    lm_fit <- lm(annulus_enrichment ~ cluster_size, data=d_sub)
    seg_fit <- segmented(lm_fit, npsi=1)
    xs <- seq(0, max(d_sub$cluster_size), length.out=1000)
    plot(d_sub$annulus_enrichment ~ d_sub$cluster_size, pch=19, xlab="Cluster size", ylab="Mean enrichment", main=tmp_subject)
    lines(predict(seg_fit, newdata=data.frame(cluster_size=xs)) ~ xs, col="red")
}
dev.off()

if(!exists("fg_dat")) {
    source("R/load_score_datasets.R")
}


p2 <- motif_dat %>%
    ggplot(aes(x=cluster_size, y=mean_enrichment, color=radius)) +
    geom_point() +
    facet_wrap(vars(subject), scales="free") +
    expand_limits(y=0) +
    scale_color_viridis()

subjects <- motif_dat[["subject"]] %>% unique
radii <- motif_dat[["radius"]] %>% unique
cluster_tcrs <- {}

if(!exists("json_object")) {
    json_object <- fromJSON(file=file.path(json_dir, "motif.json"))
}



subjects <- subjects %>% sapply(gsub, pattern=".tcrs", replacement="")
motif_rates <- list()
reverse_motif_rates <- list()
reverse_all_motif_rates <- list()
all_motif_rates <- list()
joint_rates <- list()
all_joint_rates <- list()
for(tmp_subject in subjects) {
    snap_dat <- fg_dat %>%
        build_mds_dataframe(radius=50.5, subjects=tmp_subject)
    snap_dat[["is_in_cluster"]] <- snap_dat[["tcr"]] %>% 
        sapply(function(x) { x %in% cluster_tcrs[[paste(tmp_subject, "tcrs", sep=".")]] })
    snap_dat[["has_hmmer_motif"]] <- snap_dat[["tcr"]] %>% 
        sapply(toString) %>% 
        unname %>%
        sapply(has_ida_plus_motif, e_value_dat=e_value_dat, tmp_subject=tmp_subject, e_value_threshold=1e-7)
    snap_dat %>%
        ggplot(aes(x=x1, y=x2, color=score)) + 
        geom_point(aes(shape=is_in_cluster)) +
        scale_color_viridis_c() +
        theme_minimal()
    ggsave(file.path(motif_dir, paste0(tmp_subject, ".pdf")), width=8, height=8)

    snap_dat %>%
        ggplot(aes(x=x1, y=x2, color=score)) + 
        geom_point(aes(shape=has_hmmer_motif)) +
        scale_color_viridis_c() + 
        theme_minimal()
    ggsave(file.path(e_value_motif_dir, paste0(tmp_subject, ".pdf")), width=8, height=8)

    motif_rates[[tmp_subject]] <- snap_dat[snap_dat[["label"]] %in% c("Revere", "Tremont"), ][["is_in_cluster"]] %>% mean
    all_motif_rates[[tmp_subject]] <- snap_dat[snap_dat[["label"]] %in% c("Ida", "Revere", "Tremont"), ][["is_in_cluster"]] %>% mean
    reverse_motif_rates[[tmp_subject]] <- snap_dat[snap_dat[["is_in_cluster"]], ][["label"]] %>% 
        sapply(function(x) { x %in% c("Revere", "Tremont")}) %>% 
        mean
    reverse_all_motif_rates[[tmp_subject]] <- snap_dat[snap_dat[["is_in_cluster"]], ][["label"]] %>% 
        sapply(function(x) { x %in% c("Revere", "Tremont", "Ida")}) %>% 
        mean
    joint_rates[[tmp_subject]] <- Map(function(x, y) { x && y },
                                      snap_dat[["is_in_cluster"]],
                                      (snap_dat[["label"]] %in% c("Revere", "Tremont"))
                                     ) %>% unlist %>% mean
    all_joint_rates[[tmp_subject]] <- Map(function(x, y) { x && y },
                                      snap_dat[["is_in_cluster"]],
                                      (snap_dat[["label"]] %in% c("Revere", "Tremont", "Ida"))
                                     ) %>% unlist %>% mean
}

