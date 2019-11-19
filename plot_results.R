library(dplyr)
library(ggplot2)
library(reshape2)
library(rjson)

results <- fromJSON(file="results_data.json")
names(results) <- c("S1-S1 (control)",
                    "S1-S2",
                    "S1-S3",
                    "S1-S4"
                   )
scores <- results %>%
    lapply(function(x) { x[['max_scores']] }) %>%
    melt %>%
    data.table::setnames(c("Score", "Group"))
jaccards <- results %>% 
    lapply(function(x) { x[['jaccard_indices']] }) %>%
    melt %>%
    data.table::setnames(c("Jaccard", "Group"))

p1 <- ggplot(scores) + 
    geom_density(aes(x=Score, color=Group), size=1) +
    ggtitle("Max score by group") +
    theme_minimal(base_size=14)
p2 <- ggplot(jaccards) + 
    geom_density(aes(x=Jaccard, color=Group), size=1)  +
    ggtitle("Jaccard index by group") +
    theme_minimal(base_size=14)

merged <- cbind(scores, jaccards)[, 1:3]
merged[["Fit"]] <- lm(Score ~ Jaccard, data=merged)[["fitted.values"]]
p3 <- ggplot(merged, aes(x=Jaccard, y=Score, color=Group)) + 
    geom_point() +
    geom_line(aes(x=Jaccard, y=Fit), color="black", size=1) +
    ggtitle("Max score versus Jaccard index (with least-squares fit)") +
    theme_minimal(base_size=14)


