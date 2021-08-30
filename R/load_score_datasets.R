library(rjson)
library(dplyr)
library(reshape2)

# source("R/plot_utils.R")

json_dir = CONFIG$JSON_OUTPUT
dist_mats_dir = CONFIG$DIST_MATRICES_OUTPUT
projection_method <- "MDS"
dir.create(projection_method)

if(!exists("dat")) {
    results <- list("results"=fromJSON(file=file.path(json_dir, "empirical_fg_bg_nbhd_stats.json")))
    for(subject in names(results)) {
        results[[subject]][["sds"]] <- NULL
    }
    dat <- results %>% melt
}

if(!exists("fg_dat")) {

    dat[, c("L1")] <- NULL
    names(dat) <- c("Value", "ValueType", "Group", "TCR", "TCRDistRadius", "Subject")
    neighbor_counts <- dat[dat[["ValueType"]] == "neighbor_count", ][["Value"]]
    dat <- dat[dat[["ValueType"]] == "score", ]
    dat[["NeighborCount"]] <- neighbor_counts
    names(dat)[1] <- "Score"
    dat[["TCRDistRadius"]] <- dat[["TCRDistRadius"]] %>% as.numeric
    dat[["Subject"]] <- dat[["Subject"]] %>% 
        sapply(gsub, pattern=".tcrs", replacement="")
    
    fg_dat <- dat[dat$Group == "foreground", ]
    bg_dat <- dat[dat$Group == "background", ]
    
    subjects <- dat$Subject %>% 
        unique
}

dist_mats <- {}
for(subject in subjects) {
    dist_mats[[subject]] <- read.csv(
        file.path(
            dist_mats_dir,
            paste0(
                 subject,
                 ".tcrs.csv"
            )
        ),
        header=FALSE
    )
}

mds_dats <- {}
for(subject in subjects) {
    mds_dats[[subject]] <- get_projection(dist_mats[[subject]], 
                                          projection_method=projection_method, 
                                          k=2)
}
