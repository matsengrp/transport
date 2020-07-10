library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(rjson)
library(segmented)
library(viridis)

source("R/plot_utils.R")

get_breakpoint_from_model <- function(seg_fit, default_bp=50) {
    bp <- seg_fit[["psi"]][2]
    if(is.null(bp)) {
        bp <- default_bp
    }

    return(bp)
}

csv_dir <- "output/csv"

seg_dat <- fread(file.path(csv_dir, "seg.csv"))

radii <- seg_dat[["radius"]]
lm_fit <- lm(annulus_enrichment ~ radius, data=seg_dat)
seg_fit <- segmented(lm_fit, npsi=1)
breakpoint <- seg_fit %>% get_breakpoint_from_model
cutoff_radius <- radii[which(radii <= breakpoint) %>% max]

write(cutoff_radius, file="tmp_output/breakpoint.txt")

args = commandArgs(trailingOnly=TRUE)
xs <- seq(0, max(seg_dat$radius), length.out=1000)
pdf(paste0('output/hmm/', args[1], '.pdf'))
plot(seg_dat$annulus_enrichment ~ seg_dat$radius)
lines(predict(seg_fit, newdata=data.frame(radius=xs)) ~ xs, col="red")
dev.off()
