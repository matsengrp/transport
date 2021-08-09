library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
source("R/plot_utils.R")
dir <- CONFIG$HMM_CD4_DN_OUTPUT

clusters <- 1:3 %>% sapply(function(x) { paste("cluster", x, sep="_") })
for(cluster in clusters) {
    cluster_dir <- file.path(dir, cluster)
    cluster_df <- fread(file.path(cluster_dir, "cluster_tcrs.csv"), header=TRUE)
    cluster_df[["v_gene"]] <- cluster_df[["tcr"]] %>%
        sapply(function(x) { s <- strsplit(x, split=",") %>% unlist %>% first }) %>%
        sapply(gsub, pattern="TRBV", replacement="") %>%
        sapply(as.factor)

    v_gene_freqs <- table(cluster_df[["v_gene"]])/nrow(cluster_df)

    colors <- rep(brewer.pal(12, "Paired"), 2)

    v_gene_freqs %>%
        data.table %>%
        ggplot(aes(fill=reorder(V1, -N), y=N, x=factor(""))) + 
        geom_bar(stat="identity", position="fill") +
        xlab("V gene") +
        ylab("Frequency") +
        theme_minimal() +
        theme(legend.title = element_blank(), text = element_text(size=20)) +
        scale_fill_manual(values=colors) +
    ggsave(file.path(cluster_dir, "v_genes.pdf"), width=3, height=6)
}
