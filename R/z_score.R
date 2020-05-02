library(cowplot)
library(dplyr)
library(ggplot2)
library(reshape2)
library(rjson)

source("plot_utils.R")

dn_subject <- 'DN_18_B.tcrs'

results <- fromJSON(file='replicate_z_scores.json')
df <- results %>% melt
names(df) <- c("z_score", "group", "tcr", "subject")
replicate_df <- df[df[['subject']] == dn_subject, ]
replicate_df[["tmp"]] <- NULL

rand_results <- fromJSON(file='rand_z_scores.json')
rand_df <- rand_results %>% melt
names(rand_df) <- c("z_score", "tcr")
rand_df[["group"]] <- 'randomization'

tall_df <- rbind.data.frame(
    replicate_df[, c('z_score', 'tcr', 'group')], 
    rand_df
)


names(rand_df)[which(names(rand_df) == 'z_score')] <- 'rand_z_score'

bg_df <- replicate_df[replicate_df[['group']] == 'background', c('z_score', 'tcr')]
names(bg_df)[which(names(bg_df) == 'z_score')] <- 'bg_z_score'
fg_df <- replicate_df[replicate_df[['group']] == 'foreground', c('z_score', 'tcr')]
names(fg_df)[which(names(fg_df) == 'z_score')] <- 'fg_z_score'


wide_df <- bg_df %>%
    merge(rand_df, by='tcr')

wide_df[['label']] <- wide_df[['tcr']] %>%
    sapply(get_motif_label)
tall_df[['label']] <- tall_df[['tcr']] %>%
    sapply(get_motif_label)

p_scatter <- wide_df %>%
    ggplot() +
    geom_point(aes(x=rand_z_score, y=bg_z_score, color=label, shape=label), alpha=0.6) +
    geom_smooth(method="lm", aes(y=bg_z_score, x=rand_z_score)) +
    theme_minimal() +
    xlab("Randomization z-score") +
    ylab("Background z-score")
ggsave("z_score_scatterplot.pdf")

p_densities <- tall_df[tall_df[['group']] != 'foreground', ] %>%
    ggplot(aes(x=z_score, y=..density.., color=label)) +
    facet_wrap(vars(group), dir="v", scales="free") +
    theme_minimal() +
    geom_density()
ggsave("z_score_densities.pdf")
