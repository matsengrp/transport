library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(rjson)
library(segmented)
library(viridis)

#source("R/plot_utils.R")

default_bp <- 75

get_breakpoint_from_model <- function(seg_fit, default_bp=default_bp) {
    tryCatch({
        bp <- seg_fit[["psi"]][2]
    }, error=function(cond) {})
    if(is.null(bp)) {
        bp <- default_bp
    }

    return(bp)
}

csv_dir <- "output/csv"

seg_dat <- fread(file.path(csv_dir, "seg.csv"))

radii <- seg_dat[["radius"]]
tryCatch({
        lm_fit <- lm(annulus_enrichment ~ radius, data=seg_dat)
        seg_fit <- segmented(lm_fit, npsi=1)
        breakpoint <- seg_fit %>% get_breakpoint_from_model
    }, 
    error=function(cond){ 
    }, 
    finally={
        breakpoint <- default_bp
    }
)
cutoff_radius <- radii[which(radii <= breakpoint) %>% max]

args <- commandArgs(trailingOnly=TRUE)
outdir <- args[1]
write(cutoff_radius, file=file.path(outdir, "breakpoint.txt"))

xs <- seq(0, max(seg_dat$radius), length.out=1000)
pdf(file.path(outdir, 'seg_reg.pdf'), width=8, height=4)
plot(seg_dat$annulus_enrichment ~ seg_dat$radius, pch=19, xlab="Radius", ylab="Mean loneliness")
lines(predict(seg_fit, newdata=data.frame(radius=xs)) ~ xs, col="red")
dev.off()
