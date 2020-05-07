library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(segmented)
library(viridis)

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
par(mfrow=c(4, 4))
for(tmp_subject in subjects) {
    d_sub <- dat[dat$subject == tmp_subject, ]
    lm_fit <- lm(mean_enrichment ~ radius, data=d_sub)
    seg_fit <- segmented(lm_fit)
    plot(d_sub$mean_enrichment ~ d_sub$radius, pch=19)
    lines(seg_fit$fitted.values ~ d_sub$radius, col="red")
}
