
if (config$threads > 1) {
    library(doMC)
    registerDoMC(cores=config$threads)
    do_parallel <- T
} else {
    do_parallel <- F
}


# Munge and cache correlation matrices
cache("avana_dep_corr", {
    avana_dep_corr <- avana_dep[avana_analysis_genes, ] %>%
        row_correlations(.parallel=do_parallel, use="pairwise")
}, depends = c("avana_dep", "avana_analysis_genes"))


cache("rnai_dep_corr", {
    keep_genes <- rnai_analysis_genes %>%
        intersect(rownames(rnai_dep)[aaply(rnai_dep, 1, function(x) sum(is.na(x))) < 285]) %>%
        unique

    rnai_dep_corr <- rnai_dep[keep_genes, ] %>%
        row_correlations(.parallel=do_parallel, use="pairwise")
}, depends = c("rnai_dep", "rnai_analysis_genes"))

cache("gecko_dep_corr", {
    gecko_dep[intersect(rownames(gecko_dep), avana_analysis_genes),] %>%
        row_correlations(.parallel=do_parallel, use="pairwise")
}, depends = c("gecko_dep", "avana_analysis_genes"))


cache("wang_dep_corr", {
    wang_dep[intersect(rownames(wang_dep), avana_analysis_genes),] %>%
        row_correlations(.parallel=do_parallel, use="pairwise")
}, depends = c("wang_dep", "avana_analysis_genes"))



# Munge and cache symmetrized rank matrices
cache("avana_dep_rank", {
    edgeweight_symmetric_rank(avana_dep_corr)
}, depends = "avana_dep_corr")

cache("rnai_dep_rank", {
    edgeweight_symmetric_rank(rnai_dep_corr)
}, depends = "rnai_dep_corr")

cache("gecko_dep_rank", {
    edgeweight_symmetric_rank(gecko_dep_corr)
}, depends = "wang_dep_corr")

cache("wang_dep_rank", {
    edgeweight_symmetric_rank(wang_dep_corr)
}, depends = "wang_dep_corr")



registerDoSEQ()
