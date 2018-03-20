# calculate_edgeweights(dat, include_genes=NULL) {
#
# }

edgeweight_empirical_fdr <- function(mat, null_genes, by_row=T, do_parallel=F) {
    null_genes <- intersect(null_genes, colnames(mat))
    aaply(mat, 1, function(x) {
        null_scores_ecdf <- ecdf(x[null_genes])
        number_of_nulls <- null_scores_ecdf(x)*length(null_genes)
        frac_null <- (length(null_genes) - number_of_nulls + 1)/ (length(null_genes) + 1)
        frac_null %>% set_names(names(x))
    }, .parallel=do_parallel) %>%
        `diag<-`(1)
}

edgeweight_symmetric_rank <- function(mat, rank_method = "maximal") {
    edge_rank <- mat
    diag(edge_rank) <- NA
    edge_rank <- matrixStats::rowRanks(-edge_rank)  %>%
        set_rownames(rownames(edge_rank)) %>%
        set_colnames(colnames(edge_rank))

    if (rank_method=="maximal") {
        edge_rank <- pmin(edge_rank, t(edge_rank), na.rm=T)
    }
    if (rank_method=="minimal") {
        edge_rank <- pmax(edge_rank, t(edge_rank), na.rm=T)
    }
    if (rank_method=="average") {
        edge_rank <- (edge_rank + t(edge_rank)) / 2
    }
    return(edge_rank)
}

filter_edgeweight_matrix <- function(edgeweight_mat, filter_mat, cutoff=0.01, new_val=0) {
    edgeweight_mat[filter_mat > cutoff] <- new_val
    return(edgeweight_mat)
}
