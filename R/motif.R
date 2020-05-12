library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(rjson)
library(segmented)
library(viridis)

get_breakpoint_from_model <- function(seg_fit) {
    return(seg_fit[["psi"]][2])
}

csv_dir <- "output/csv"
json_dir <- "output/json"

dat <- fread(file.path(csv_dir, "motif.csv"))

p1 <- dat %>%
    ggplot(aes(x=radius, y=mean_enrichment, label=cluster_size)) +
    geom_text(size=2) +
    facet_wrap(vars(subject), scales="free") +
    expand_limits(y=0) +
    scale_color_viridis()

p2 <- dat %>%
    ggplot(aes(x=cluster_size, y=mean_enrichment, color=radius)) +
    geom_point() +
    facet_wrap(vars(subject), scales="free") +
    expand_limits(y=0) +
    scale_color_viridis()

subjects <- dat[["subject"]] %>% unique
radii <- dat[["radius"]] %>% unique
cluster_tcrs <- {}

if(!exists("json_object")) {
    json_object <- fromJSON(file=file.path(json_dir, "motif.json"))
}

par(mfrow=c(4, 4))
for(tmp_subject in subjects) {
    d_sub <- dat[dat$subject == tmp_subject, ]
    lm_fit <- lm(mean_enrichment ~ radius, data=d_sub)
    seg_fit <- segmented(lm_fit)
    breakpoint <- seg_fit %>% get_breakpoint_from_model
    cutoff_radius <- radii[which(radii < breakpoint) %>% max]
    cluster_tcrs[[tmp_subject]] <- json_object[[tmp_subject]][[toString(cutoff_radius)]][["tcrs"]]
    xs <- seq(0, max(d_sub$radius), length.out=1000)
    plot(d_sub$mean_enrichment ~ d_sub$radius, pch=19, xlab="Radius", ylab="Mean enrichment", main=tmp_subject)
    lines(predict(seg_fit, newdata=data.frame(radius=xs)) ~ xs, col="red")
}

if(!exists("fg_dat")) {
    source("R/load_score_datasets.R")
}

motif_dir <- "output/motif"
dir.create(motif_dir)



subjects <- subjects %>% sapply(gsub, pattern=".tcrs", replacement="")
for(tmp_subject in subjects) {
    snap_dat <- fg_dat %>%
        build_mds_dataframe(radius=50.5, subjects=tmp_subject)
    snap_dat[["is_in_cluster"]] <- snap_dat[["tcr"]] %>% 
        sapply(toString) %>% 
        sapply(function(x) { x %in% cluster_tcrs[[paste(tmp_subject, "tcrs", sep=".")]] })
    snap_dat %>%
        ggplot(aes(x=x1, y=x2, color=score)) + 
        geom_point(aes(shape=is_in_cluster)) +
        scale_color_viridis_c()
    ggsave(file.path(motif_dir, paste0(tmp_subject, ".pdf")), width=8, height=8)
}

