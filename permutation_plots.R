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

process_results <- function(cell_type) {
	results <- fromJSON(file=file.path(dir, cell_type, "within_results.json"))
	cell_type_df <- results %>% melt(stringsAsFactors=TRUE)
	cell_type_df <- cell_type_df[, c(1, 2, 4, 5, 6, 7, 8, 10)]
	names(cell_type_df) <- c("Score", "CDR3", "Margin", "Gene", "Subject", "Trial", "DistributionType", "Group")
	return(cell_type_df)
}

merge_dfs <- function(cell_df, permutation_df, cell_type, margin="column") {
	group_to_remove_hash <- list("CD4"="CD8_DN",
								 "CD8"="CD4_DN",
								 "DN"="CD4_CD8")
	columns_to_keep <- c("Value", "ValueType", "Margin", "Gene", "Group", "DistributionType", "Trial")
	full_df <- rbind.data.frame(cell_df[, columns_to_keep], 
								permutation_df[, columns_to_keep])
	df_sub <- full_df[full_df$Margin == margin & 
					 full_df$ValueType == "scores" &
	                 full_df$Group != group_to_remove_hash[[cell_type]], ]
	if(cell_type == "CD4") {
		df_sub[df_sub$Group == "CD4_CD8", ]$Group <- "CD8"
		df_sub[df_sub$Group == "CD4_DN", ]$Group <- "DN"
	} else if(cell_type == "CD8") {
		df_sub[df_sub$Group == "CD4_CD8", ]$Group <- "CD4"
		df_sub[df_sub$Group == "CD8_DN", ]$Group <- "DN"
	} else if(cell_type == "DN") {
		df_sub[df_sub$Group == "CD4_DN", ]$Group <- "CD4"
		df_sub[df_sub$Group == "CD8_DN", ]$Group <- "CD8"
	}

	return(df_sub)
}

make_cell_plots <- function(cell_df, permutation_df, cell_type) {
	dir.create(file.path(dir, cell_type))

	df <- merge_dfs(cell_df, permutation_df, cell_type)
	cell_plot <- (df %>% 
		ggplot(aes(x=Value, colour=Group)) +
		ggtitle(paste("Aggregate CDFs for reference group =", cell_type)) +
		stat_ecdf(aes(lty=DistributionType))) %>%
		format_plot
	ggsave(file.path(dir, cell_type, paste0("all_ecdfs_", tolower(cell_type), ".pdf")))

	per_trial_cell_plot <- (df %>% 
		ggplot() + 
		ggtitle(paste("Per-trial CDFs for reference group =", cell_type)) +
		stat_ecdf(aes(x=Value, group=interaction(Trial, Group), colour=Group, linetype=DistributionType), size=0.1) ) %>% 
		format_plot
	ggsave(file.path(dir, cell_type, "per_trial_permutaiton_ecdfs.pdf"))

	per_gene_cell_plot <- 
		(df %>% 
			ggplot() + 
			ggtitle(paste("Per-gene per-trial CDFs for reference group =", cell_type)) +
			stat_ecdf(aes(x=Value, group=interaction(Trial, Group), colour=Group, linetype=DistributionType), size=0.1) +
			facet_wrap(vars(Gene))
		) %>% format_plot
	ggsave(file.path(dir, cell_type, paste0("per_genes_", cell_type, ".pdf")))

	per_gene_cell_plot <- 
		(df %>% ggplot(aes(x=Value, colour=Group)) + 
				stat_ecdf(aes(lty=DistributionType), size=0.1) +
				ggtitle(paste("Aggregate per-gene CDFs for reference group =", cell_type)) +
				facet_wrap(vars(Gene))
		) %>% format_plot
	ggsave(file.path(dir, cell_type, paste0("aggregate_per_genes_", tolower(cell_type), ".pdf")))
}
	

method <- "subsample"
dir <- file.path("~/sync", method)
cd4_df <- "CD4" %>% process_results
cd8_df <- "CD8" %>% process_results
dn_df <- "DN" %>% process_results
stop()
permutation_results <- fromJSON(file=file.path(dir, "permutation_results.json"))
permutation_df <- permutation_results %>% melt(stringsAsFactors=TRUE)
permutation_df <- permutation_df[, c(1, 2, 4, 5, 6, 7, 9)]
names(permutation_df) <- c("Score", "CDR3", "Margin", "Gene", "Trial", "DistributionType", "Group")
permutation_df[["Trial"]] <- as.factor(permutation_df[["Trial"]])
permutation_df[["Group"]] <- as.factor(permutation_df[["Group"]])


make_cell_plots(cd4_df, permutation_df, "CD4")
make_cell_plots(cd8_df, permutation_df, "CD8")
make_cell_plots(dn_df, permutation_df, "DN")
