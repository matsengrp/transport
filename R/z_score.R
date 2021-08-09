library(cowplot)
library(data.table)
library(dplyr)
library(estimatr)
library(ggplot2)
library(RcmdrMisc)
library(reshape2)
library(rjson)

source("R/plot_utils.R")

lm_eqn <- function(df){
    m <- summarySandwich(lm(bg_z_score ~ rand_z_score, df));
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
             r2 = format(m$r.squared, digits = 3)))
    as.character(as.expression(eq));
}


json_dir = CONFIG$JSON_OUTPUT
z_score_dir = CONFIG$ZSCORE_OUTPUT

dn_subject <- 'DN_15_B'
cd4_subject <- 'CD4_17_B'

motif_levels <- c("N/A", "Ida", "Revere", "Tremont")
cluster_names <- list(
                      "0"="N/A",
                      "1"="OT-Tremont", 
                      "2"="OT-Revere",
                      "3"="OT-Ida",
                      "4"="N/A",
                      "5"="N/A"
                     )

cluster_df <- fread(file.path(json_dir, gsub(cd4_subject, pattern="_B", replacement=""), "cluster_df.csv"))
cluster_df[["motif"]] <- factor(cluster_df[["motif"]], levels=motif_levels)
cluster_df[["cluster"]] <- cluster_df[["cluster"]] %>%
    sapply(function(x) { cluster_names[[toString(x)]] })

results <- fromJSON(file=file.path(json_dir, 'replicate_z_scores.json'))
df <- results %>% melt
names(df) <- c("z_score", "group", "tcr", "subject")
replicate_df <- df[df[['subject']] == paste0(dn_subject, ".tcrs"), ]
replicate_df[["tmp"]] <- NULL
replicate_df[["cluster"]] <- rep(cluster_df[["cluster"]], 1, each=2) %>% sapply(toString)
replicate_df[["motif"]] <- rep(cluster_df[["motif"]], 1, each=2) %>% sapply(toString) %>%
    factor(levels=motif_levels)

rand_results <- fromJSON(file=file.path(json_dir, gsub(cd4_subject, pattern="_B", replacement=""), 'rand_z_scores.json'))
rand_df <- rand_results %>% melt
names(rand_df) <- c("z_score", "tcr")
rand_df[["group"]] <- 'randomization'
rand_df[["cluster"]] <- cluster_df[["cluster"]] %>% sapply(toString)
rand_df[["motif"]] <- cluster_df[["motif"]] %>% sapply(toString) %>%
    factor(levels=motif_levels)


tall_df <- rbind.data.frame(
    replicate_df[, c('z_score', 'tcr', 'group', 'cluster', 'motif')], 
    rand_df
)

names(rand_df)[which(names(rand_df) == 'z_score')] <- 'rand_z_score'

bg_df <- replicate_df[replicate_df[['group']] == 'background', c('z_score', 'tcr')]
names(bg_df)[which(names(bg_df) == 'z_score')] <- 'bg_z_score'
fg_df <- replicate_df[replicate_df[['group']] == 'foreground', c('z_score', 'tcr')]
names(fg_df)[which(names(fg_df) == 'z_score')] <- 'fg_z_score'

wide_df <- bg_df %>%
    merge(rand_df, by='tcr')

p_scatter <- wide_df %>%
    ggplot() +
    geom_point(aes(x=rand_z_score, y=bg_z_score, color=cluster, shape=cluster), alpha=0.6) +
    geom_smooth(method="lm_robust", aes(y=bg_z_score, x=rand_z_score)) +
    theme_minimal() +
    geom_text(x=-1, y=8, label=lm_eqn(wide_df), parse = TRUE) +
    xlab("Randomization z-score") +
    ylab("Replicate z-score")
ggsave(file.path(z_score_dir, "z_score_scatterplot.pdf"), width=8, height=6)


tall_df[tall_df[["group"]] == "background", ][["group"]] <- "replicate"
p_densities <- tall_df[tall_df[['group']] != 'foreground', ] %>%
    ggplot(aes(x=z_score, y=..density.., color=cluster)) +
    facet_wrap(vars(group), dir="v") +
    xlim(-5, 12) +
    theme_minimal() +
    geom_density()
ggsave(file.path(z_score_dir, "z_score_densities.pdf"), width=8, height=6)
