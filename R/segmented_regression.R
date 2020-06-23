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

seg_dat <- fread(file.path(csv_dir, "seg.csv"))

radii <- seg_dat[["radius"]]
lm_fit <- lm(annulus_enrichment ~ radius, data=seg_dat)
seg_fit <- segmented(lm_fit, psi=list(radius=20.5))
breakpoint <- seg_fit %>% get_breakpoint_from_model
cutoff_radius <- radii[which(radii <= breakpoint) %>% max]

write(cutoff_radius, file="tmp_output/breakpoint.txt")
