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


has_ida_plus_motif <- function(tcr, e_value_dat, tmp_subject, e_value_threshold=1e-5) {
    hmmer_cdr3s <- e_value_dat[e_value_dat[["subject"]] == paste0(tmp_subject, ".tcrs") & e_value_dat[["e_value"]] < e_value_threshold, ][["cdr3"]]
    tcr_info <- tcr %>%
         extract_tcr_info
    cdr3 <- tcr_info[2]
    return(cdr3 %in% hmmer_cdr3s)
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

get_motif_label <- function(tcr) {
    if(tcr %>% has_revere_motif) {
        label <- "Revere"
    } else if(tcr %>% has_tremont_motif) {
        label <- "Tremont"
    } else if(tcr %>% has_ida_motif) {
        label <- "Ida"
    } else if(tcr %>% has_motif_x) {
        label <- "X"
    } else {
        label <- "N/A"
    }
    return(label)
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

extract_tcrs_from_mds_cluster <- function(
                                          ref_dat,
                                          subject,
                                          xrange,
                                          yrange,
                                          score_threshold=1800,
                                          radius=50.5,
                                          score_column="score"
                                         ) {
   tmp <- ref_dat %>%
       build_mds_dataframe(radius=radius, subjects=subject)
   filtered_dat <- tmp[tmp[["x1"]] > xrange[1] &
                       tmp[["x1"]] < xrange[2] &
                       tmp[["x2"]] > yrange[1] &
                       tmp[["x2"]] < yrange[2] &
                       tmp[[score_column]] > score_threshold,
                      ]
   tcrs <- dat[dat[["Subject"]] == subject, ][["TCR"]] %>% 
       rle %$%
       values
   cluster_tcrs <- tcrs[filtered_dat %>% rownames %>% as.numeric]
   return(cluster_tcrs)
}

build_mds_dataframe <- function(ref_dat, radius, subjects, add_extra_metrics=FALSE, e_value_dat) {
    full_dat <- matrix(NA, nrow=0, ncol=4) %>% 
        data.frame %>%
        setNames(c("x1", "x2", "subject", "score"))
    for(subject in subjects) {
        subject_radius_dat <- ref_dat[ref_dat[["TCRDistRadius"]] == radius & 
                                      ref_dat[["Subject"]] == subject, ]
        subject_dat <- data.frame(
                                  x1=mds_dats[[subject]][, 1], 
                                  x2=mds_dats[[subject]][, 2],
                                  subject=subject,
                                  score=subject_radius_dat[["Score"]],
                                  label=subject_radius_dat[["TCR"]] %>% sapply(get_motif_label),
                                  tcr=subject_radius_dat[["TCR"]],
                                  revere=subject_radius_dat[["TCR"]] %>% sapply(has_revere_motif),
                                  tremont=subject_radius_dat[["TCR"]] %>% sapply(has_tremont_motif),
                                  ida=subject_radius_dat[["TCR"]] %>% sapply(has_ida_motif)
                                 )
        if(!missing(e_value_dat)) {
             subject_dat[["ida_plus"]] <- subject_radius_dat[["TCR"]] %>% 
                 sapply(
                     has_ida_plus_motif,
                     e_value_dat=e_value_dat,
                     tmp_subject=subject
                 )
             subject_dat[["ida_plus_plus"]] <- subject_radius_dat[["TCR"]] %>% 
                 sapply(
                     has_ida_plus_motif,
                     e_value_dat=e_value_dat,
                     tmp_subject=subject,
                     e_value_threshold=1e-8
                 )

        }
        subject_dat[["label"]] <- factor(subject_dat[["label"]], levels=c("N/A", "Revere", "Tremont", "Ida", "X"))

        if(add_extra_metrics) {
            subject_dat[["relative_score"]] <- subject_dat[["score"]]/max(subject_dat[["score"]])
            label_freqs <- subject_dat[["label"]] %>% 
                table %>% 
                sapply(function(x) { x/nrow(subject_dat) })
            subject_dat[["prevalence"]] <- subject_dat[["label"]] %>%
                sapply(function(x) { label_freqs[x] })
    
            ecdf_function <- ecdf(subject_dat[["score"]])
            subject_dat[["ecdf"]] <- subject_dat[["score"]] %>%
                ecdf_function
        }
        full_dat <- rbind.data.frame(full_dat, subject_dat) 
    }
    return(full_dat)
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

