library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(rjson)
library(segmented)
library(viridis)

library("rjson")
CONFIG <- fromJSON(file = "config.json")

get_breakpoint_from_model <- function(seg_fit) {
    return(seg_fit[["psi"]][2])
}

csv_dir <- CONFIG$CSV_OUTPUT
json_dir <- CONFIG$JSON_OUTPUT

motif_dir <- CONFIG$MOTIF_OUTPUT
motif_dat <- fread(file.path(csv_dir, "motif.csv"))
score_results <- list("results"=fromJSON(file=file.path(json_dir, "empirical_fg_bg_nbhd_stats.json"))) 
subjects <- motif_dat[["subject"]] %>% unique
radii <- motif_dat[["radius"]] %>% unique


breakpoints <- {}
pdf(file.path(motif_dir, "enrichment_by_radius.pdf"), width=10, height=12)
par(mfrow=c(6, 4))
motif_dat[["subject"]] <- motif_dat[["subject"]] %>% sapply(gsub, pattern=".tcrs", replacement="")
# Enforce an ordering for the subjects for plotting
subjects <- factor(subjects, levels=1:23) %>% paste("DN", ., "B", sep="_")
#print(subjects)
for(tmp_subject in subjects[order(subjects)]) {
    d_sub <- motif_dat[motif_dat$subject == tmp_subject, ]
    sample_size <- length(score_results[['results']][[paste0(tmp_subject, '.tcrs')]][['48.5']])
    #print(sample_size)
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
