dat_sub <- dat[, c("Group", "Subject", "Gene", "Count")]
dat_sub <- dat_sub[!duplicated(dat_sub), ]

print(paste("CD4:", dat_sub[dat_sub$Group == "CD4", ] %>% nrow))
print(paste("CD8:", dat_sub[dat_sub$Group == "CD8", ] %>% nrow))
print(paste("DN:", dat_sub[dat_sub$Group == "DN", ] %>% nrow))
