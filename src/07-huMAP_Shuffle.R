library(igraph)
library(ggridges)

library(ProjectTemplate)
load.project(override.config = list(munging=F, cache_loading=F))


if (config$threads > 1) {
    library(doMC)
    registerDoMC(cores=config$threads)
    do_parallel <- T
} else {
    do_parallel <- F
}


load("./cache/humap_list.RData")
# load("./cache/avana_dep_corr.RData")
load("./cache/avana_2017_dep_corr.RData")
# load("./cache/rnai_dep_corr.RData")


# humap_avana_dat <- prepare_dataset(avana_dep_corr, humap_list)
humap_avana_2017_dat <- prepare_dataset(avana_2017_dep_corr, humap_list)
# humap_rnai_dat <- prepare_dataset(rnai_dep_corr, humap_list)



ntrials <- 10000
niter <- 10
epoch <- 1000

rank_thresholds <- 2^(6:0)

humap_shuffle_wrapper <- function(shuffle_file, corr_dat) {

    if (!file.exists(shuffle_file)) {
        cat("Running network shuffle and saving to", shuffle_file, "\n")
        shuffle_result <-
            run_int_ext_shuffle(corr_dat$edgelist,
                                corr_dat$genesets,
                                rank_thresholds,
                                ntrials=ntrials, niter=niter, epoch=epoch,
                                do_parallel=do_parallel)

        saveRDS(shuffle_result, shuffle_file)
    } else {
        cat("Found", shuffle_file, "... Not re-running\n")
        shuffle_result <- readRDS(shuffle_file)
    }
    return(shuffle_result)
}

# humap_avana_shuffle <-
#     humap_shuffle_wrapper("./data/interim/shuffle_humap_avana_10k.rds",
#                           humap_avana_dat)

humap_avana_2017_shuffle <-
    humap_shuffle_wrapper("./data/interim/shuffle_humap_avana_2017_10k.rds",
                          humap_avana_2017_dat)

# humap_rnai_shuffle <-
#     humap_shuffle_wrapper("./data/interim/shuffle_humap_rnai_10k.rds",
#                           humap_rnai_dat)
