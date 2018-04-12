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


load("./cache/corum_list.RData")
load("./cache/avana_dep_corr.RData")
load("./cache/avana_2017_dep_corr.RData")
load("./cache/gecko_dep_corr.RData")
load("./cache/coxpres_db.RData")
load("./cache/rnai_dep_corr.RData")
load("./cache/wang_dep_corr.RData")

corum_avana_dat <- prepare_dataset(avana_dep_corr, corum_list)
corum_avana_2017_dat <- prepare_dataset(avana_2017_dep_corr, corum_list)
corum_rnai_dat <- prepare_dataset(rnai_dep_corr, corum_list)
corum_gecko_dat <- prepare_dataset(gecko_corr, corum_list)
corum_wang_dat <- prepare_dataset(wang_corr, corum_list)
corum_coxpr_dat <- prepare_dataset(coxpres_db, corum_list)

ntrials <- 10000
niter <- 10
epoch <- 1000

rank_thresholds <- 2^(8:0)

corum_shuffle_wrapper <- function(shuffle_file, corr_dat) {
    if (!file.exists(shuffle_file)) {
        shuffle_result <-
            run_int_ext_shuffle(corr_dat$edgelist,
                                corr_dat$genesets,
                                rank_thresholds,
                                ntrials=ntrials, niter=niter, epoch=epoch,
                                do_parallel=do_parallel)

        saveRDS(shuffle_result,
                shuffle_file)
    } else {
        shuffle_result <- readRDS(shuffle_file)
    }
    return(shuffle_result)
}


corum_avana_shuffle <-
    corum_shuffle_wrapper("./data/interim/shuffle_corum_avana_10k.rds",
                          corum_avana_dat)

corum_avana_2017_shuffle <-
    corum_shuffle_wrapper("./data/interim/shuffle_corum_avana_2017_10k.rds",
                          corum_avana_2017_dat)

corum_rnai_shuffle <-
    corum_shuffle_wrapper("./data/interim/shuffle_corum_rnai_10k.rds",
                          rnai_dep_corr)

corum_gecko_shuffle <-
    corum_shuffle_wrapper("./data/interim/shuffle_corum_gecko_10k.rds",
                          gecko_dep_corr)

corum_wang_shuffle <-
    corum_shuffle_wrapper("./data/interim/shuffle_corum_wang_10k.rds",
                          wang_dep_corr)

corum_coxpr_shuffle <-
    corum_shuffle_wrapper("./data/interim/shuffle_corum_coxpr_10k.rds",
                          coxpres_db)

