#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(cowplot)
library(dplyr)
library(ggplot2)
library(reshape2)
library(rjson)

format_plot <- function(p) {
    formatted_p <- p +
        theme_minimal(base_size=14) +
        theme(axis.text.x = element_text(size=14),
              axis.text.y = element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_text(size=16),
              plot.title = element_text(hjust = 0.5)
             )
    return(formatted_p)    
}

method <- "subsample"
reference_group <- args[1] #"CD8"
results <- fromJSON(file=file.path(method, reference_group, "within_results.json"))
names(results) <- c(
                    "CD4",
                    "CD8",
                    "DN"
                   )



dat = matrix(nrow=0, ncol=5) %>% data.frame %>% setNames(c("Group", "Subject", "Gene", "Margin", "Score"))
for(group_name in results %>% names) {
    group = results[[group_name]]
    for(subject_name in group %>% names) {
        print(subject_name)
        subject = group[[subject_name]]
        for(gene_name in subject %>% names) {
            gene = subject[[gene_name]]
            for(margin_name in gene %>% names) {
                score_hash = gene[[margin_name]]
                sub_dat = data.frame(Group=group_name, Subject=subject_name, Gene=gene_name, Margin=margin_name, Score=score_hash[["scores"]], Count=score_hash[["count"]])
                dat = rbind.data.frame(dat, sub_dat)
            }
        }
    }
}

dat_row_scores = dat[dat[["Margin"]] == "row", ]
dat_column_scores = dat[dat[["Margin"]] == "column", ]
p1 <- (ggplot(dat_row_scores) + 
       stat_ecdf(aes(x=Score, group=Subject, colour=Group), size=0.1) +
       ggtitle("Row score ECDFs by cell group and subject")
   ) %>% format_plot
p2 <- (ggplot(dat_row_scores) + stat_ecdf(aes(x=Score, group=Group, colour=Group))
        + ggtitle("Row score ECDFs by cell group")
   ) %>% format_plot
p3 <- (ggplot(dat_column_scores) + 
       stat_ecdf(aes(x=Score, group=Subject, colour=Group), size=0.1) +
       ggtitle("Column score ECDFs by cell group and subject")
   ) %>% format_plot
p4 <- (ggplot(dat_column_scores) + 
       stat_ecdf(aes(x=Score, group=Group, colour=Group)) +
       ggtitle("Column score ECDFs by cell group")
   ) %>% format_plot
grid_plot <- cowplot::plot_grid(p1, p2, p3, p4)
plot_dir <- file.path("~/sync", method)
plot_dir %>% dir.create
ggsave(file.path(plot_dir, paste0("ecdfs_", tolower(reference_group), ".pdf")), height=10, width=16)

p_by_gene <- (ggplot(dat[dat[["Margin"]] == "column", ]) + stat_ecdf(aes(x=Score, group=Group, colour=Group)) + facet_wrap(vars(Gene))) %>% format_plot
ggsave(file.path(plot_dir, paste0("ecdfs_", tolower(reference_group), "_by_gene.pdf")), height=10, width=16)

