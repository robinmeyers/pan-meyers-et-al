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


out_dir <- file.path("./output/empirical_shuffle", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)

load("./cache/humap_list.RData")
# load("./cache/avana_dep_corr.RData")
load("./cache/avana_2017_dep_corr.RData")
# load("./cache/rnai_dep_corr.RData")


# humap_avana_shuffle <- readRDS("./data/interim/shuffle_humap_avana_10k.rds")
humap_avana_2017_shuffle <- readRDS("./data/interim/shuffle_humap_avana_2017_10k.rds")
# humap_rnai_shuffle <- readRDS("./data/interim/shuffle_humap_rnai_10k.rds")


# humap_avana_dat <- prepare_dataset(avana_dep_corr, humap_list)
humap_avana_2017_dat <- prepare_dataset(avana_2017_dep_corr, humap_list)
# humap_rnai_dat <- prepare_dataset(rnai_dep_corr, humap_list)

rank_thresholds <- 2^(0:13)

# humap_avana_true <-
#     run_int_ext_true(humap_avana_dat$edgelist,
#                      humap_avana_dat$genesets,
#                      rank_thresholds)
#
# humap_avana_true %>%
#     edgedens_true_gg(title="hu.MAP mapped onto Avana Similarity Network") +
#     ggsave(file.path(out_dir, "edgedens_true_humap_avana.pdf"),
#            width=6, height=4)


humap_avana_2017_true <-
    run_int_ext_true(humap_avana_2017_dat$edgelist,
                     humap_avana_2017_dat$genesets,
                     rank_thresholds)

humap_avana_2017_true %>%
    edgedens_true_gg(title="hu.MAP mapped onto Avana2017 Similarity Network") +
    ggsave(file.path(out_dir, "edgedens_true_humap_avana_2017.pdf"),
           width=6, height=4)


# humap_rnai_true <-
#     run_int_ext_true(humap_rnai_dat$edgelist,
#                      humap_rnai_dat$genesets,
#                      rank_thresholds)
#
# humap_rnai_true %>%
#     edgedens_true_gg(title="hu.MAP mapped onto RNAi Similarity Network") +
#     ggsave(file.path(out_dir, "edgedens_true_humap_rnai.pdf"),
#            width=6, height=4)




humap_true <- bind_rows(
    # humap_avana_true %>% mutate(Network="Avana"),
    humap_avana_2017_true %>% mutate(Network="Avana 2017")
    # humap_rnai_true %>% mutate(Network="RNAi"),
)


humap_shuffle <- bind_rows(
    # humap_avana_shuffle %>% mutate(Network="Avana"),
    humap_avana_2017_shuffle %>% mutate(Network="Avana 2017")
    # humap_rnai_shuffle %>% mutate(Network="RNAi"),
)



empirical_shuffle_results <-
    inner_join(humap_true, humap_shuffle) %>%
    group_by(Network, Rank, Geneset) %>%
    mutate(IntExtTrue = ifelse(is.nan(IntExtTrue), 0, IntExtTrue)) %>%
    summarise(N_Above = sum(ifelse(is.nan(unlist(IntExtShuffle)),
                                   0 >= IntExtTrue,
                                   unlist(IntExtShuffle) >= IntExtTrue)),
              Total = length(unlist(IntExtShuffle))) %>%
    mutate(EmpiricalP = (N_Above+1)/(Total+1)) %>%
    group_by(Network) %>%
    mutate(FDR = p.adjust(EmpiricalP, method='fdr'))

rank_fdrs <- empirical_shuffle_results %>%
    group_by(Network, Geneset) %>%
    summarise(BestFDR = min(FDR, na.rm=T),
              BestRank = Rank[FDR == BestFDR][1],
              FirstRank = min(Rank[FDR < 0.05], na.rm=T)) %>%
    group_by(Network) %>%
    arrange(FirstRank) %>%
    mutate(Recall = cumsum(BestFDR < 0.05)/length(humap_list))

write_tsv(rank_fdrs, file.path(out_dir, "humap_shuffle_10k_results.tsv"))

rank_fdrs %>%
    filter(BestFDR < 0.05) %>%
    select(Network, Geneset) %>%
    write_tsv(file.path("data/interim/humap_shuffle_10k_sig_genesets.tsv"))



ggplot(rank_fdrs, aes(FirstRank, Recall, color=Network)) +
    geom_hline(aes(yintercept=RecallFinal, color=Network), linetype=2,
               data=rank_fdrs %>% group_by(Network) %>%
                   summarize(RecallFinal = sum(BestFDR<0.05) /
                                 length(humap_list))) +
    geom_line() +
    scale_y_continuous(limits=c(0, 0.5), expand=c(0, 0)) +
    scale_x_continuous(trans="log2", breaks=rank_thresholds) +
    labs(x="Rank Threshold", y="Recall of Complexes",
         title="hu.MAP mapped onto Similarity Networks") +
    ggsave(file.path(out_dir, "rank_recall_humap.pdf"),
           width=6, height=4)

rank_fdrs %>%
    group_by(Network) %>%
    arrange(BestFDR) %>%
    mutate(Recall = cumsum(BestFDR < 1)/length(humap_list)) %>%
    ggplot(aes(Recall, 1-BestFDR, color=Network)) +
    geom_hline(yintercept=0.95, linetype=2) +
    geom_vline(aes(xintercept=RecallFinal, color=Network), linetype=2,
               data=rank_fdrs %>% group_by(Network) %>%
                   summarize(RecallFinal = sum(BestFDR<0.05) /
                                 length(humap_list))) +
    geom_line() +
    scale_y_continuous(limits=c(0.5, 1.02), expand=c(0, 0)) +
    scale_x_continuous(expand=c(0, 0)) +
    labs(y="Precision", x="Recall",
         title="hu.MAP Precision-Recall in Similarity Networks") +
    ggsave(file.path(out_dir, "precision_recall_humap.pdf"),
           width=6, height=4)

