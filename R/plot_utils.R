
library("rjson")
CONFIG <- fromJSON(file = "config.json")
csv_dir <- CONFIG["CSV_OUTPUT"]

extract_tcr_info <- function(tcr) {
    return(tcr %>% strsplit(split=",") %>% first)
}

has_revere_motif <- function(tcr) {
    tcr_info <- tcr %>% 
        extract_tcr_info
    gene <- tcr_info[1]
    cdr3 <- tcr_info[2]
    revere_cdr3 <- "GT[VI]SNERLFF"
    return(grepl(revere_cdr3, cdr3))
}

has_tremont_motif <- function(tcr) {
    tcr_info <- tcr %>%
        extract_tcr_info
    gene <- tcr_info[1]
    cdr3 <- tcr_info[2]
    tremont_gene <- "TRBV16"
    tremont_cdr3 <- "DWG"
    return(grepl(tremont_gene, gene) && grepl(tremont_cdr3, cdr3))
}

has_ida_motif <- function(tcr) {
    tcr_info <- tcr %>%
        extract_tcr_info
    gene <- tcr_info[1]
    cdr3 <- tcr_info[2]
    ida_gene <- "TRBV(16|12)"
    ida_cdr3 <- "Y*A*EQ[YF]F"
    return(grepl(ida_gene, gene) && grepl(ida_cdr3, cdr3) && !has_tremont_motif(tcr))
}


has_motif_x <- function(tcr) {
    tcr_info <- tcr %>%
        extract_tcr_info
    gene <- tcr_info[1]
    cdr3 <- tcr_info[2]
    x_gene <- "TRBV(12|13)"
    x_cdr3 <- "[AEN]E[RT]L[FY]F"
    #x_cdr3  <- "[AEN]E[TR][LQ][FY]F"
    return(grepl(x_gene, gene) && grepl(x_cdr3, cdr3) && !has_revere_motif(tcr))
}

get_projection <- function(dist_mat, projection_method="MDS", k=2) {
    if(projection_method == "MDS") {
        projection <- cmdscale(dist_mat, k=k)
    } else if(projection_method == "TSNE") {
        library(tsne)
        projection <- tsne(dist_mat, k=k)
    } else {
        stop(paste("Unsupported projection method:", projection_method))
    } 

    return(projection)
}


plot_mds_with_scores <- function(ref_dat, radius, subject) {
    mds_dat <- build_mds_dataframe(ref_dat, radius, subject)
    p <- ggplot(mds_dat, aes(x=x1, y=x2, color=score)) + 
        geom_point(size=0.5) +
        scale_colour_viridis_c(breaks = as.numeric(mds_dat$score)) +
        theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(),
              panel.background = element_blank()) +
        ggtitle(paste("Radius = ", radius))
    return(p)
}

