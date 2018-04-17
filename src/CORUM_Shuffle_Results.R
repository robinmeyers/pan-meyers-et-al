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

load("./cache/corum_list.RData")

load("./cache/avana_dep_corr.RData")
load("./cache/avana_2017_dep_corr.RData")
load("./cache/gecko_dep_corr.RData")
load("./cache/coxpres_db.RData")
load("./cache/rnai_dep_corr.RData")
load("./cache/wang_dep_corr.RData")


corum_avana_shuffle <- readRDS("./data/interim/shuffle_corum_avana_10k.rds")
corum_avana_2017_shuffle <- readRDS("./data/interim/shuffle_corum_avana_2017_10k.rds")
corum_rnai_shuffle <- readRDS("./data/interim/shuffle_corum_rnai_10k.rds")
corum_gecko_shuffle <- readRDS("./data/interim/shuffle_corum_gecko_10k.rds")
corum_wang_shuffle <- readRDS("./data/interim/shuffle_corum_wang_10k.rds")
corum_coxpr_shuffle <- readRDS("./data/interim/shuffle_corum_coxpr_10k.rds")




corum_avana_dat <- prepare_dataset(avana_dep_corr, corum_list)
corum_avana_2017_dat <- prepare_dataset(avana_2017_dep_corr, corum_list)
corum_rnai_dat <- prepare_dataset(rnai_dep_corr, corum_list)
corum_gecko_dat <- prepare_dataset(gecko_dep_corr, corum_list)
corum_wang_dat <- prepare_dataset(wang_dep_corr, corum_list)
corum_coxpr_dat <- prepare_dataset(coxpres_db, corum_list)

rank_thresholds <- 2^(0:13)

corum_avana_true <-
    run_int_ext_true(corum_avana_dat$edgelist,
                     corum_avana_dat$genesets,
                     rank_thresholds)

corum_avana_true %>%
    edgedens_true_gg(title="CORUM mapped onto Avana Similarity Network") +
    ggsave(file.path(out_dir, "edgedens_true_corum_avana.pdf"),
           width=6, height=4)


corum_avana_2017_true <-
    run_int_ext_true(corum_avana_2017_dat$edgelist,
                     corum_avana_2017_dat$genesets,
                     rank_thresholds)

corum_avana_2017_true %>%
    edgedens_true_gg(title="CORUM mapped onto Avana2017 Similarity Network") +
    ggsave(file.path(out_dir, "edgedens_true_corum_avana_2017.pdf"),
           width=6, height=4)


corum_rnai_true <-
    run_int_ext_true(corum_rnai_dat$edgelist,
                     corum_rnai_dat$genesets,
                     rank_thresholds)

corum_rnai_true %>%
    edgedens_true_gg(title="CORUM mapped onto RNAi Similarity Network") +
    ggsave(file.path(out_dir, "edgedens_true_corum_rnai.pdf"),
           width=6, height=4)


corum_coxpr_true <-
    run_int_ext_true(corum_coxpr_dat$edgelist,
                     corum_coxpr_dat$genesets,
                     rank_thresholds)

corum_coxpr_true %>%
    edgedens_true_gg(title="CORUM mapped onto COXPRESdb Similarity Network") +
    ggsave(file.path(out_dir, "edgedens_true_corum_coexpressdb.pdf"),
           width=6, height=4)

corum_gecko_true <-
    run_int_ext_true(corum_gecko_dat$edgelist,
                     corum_gecko_dat$genesets,
                     rank_thresholds)

corum_gecko_true %>%
    edgedens_true_gg(title="CORUM mapped onto GeCKOv2 Similarity Network") +
    ggsave(file.path(out_dir, "edgedens_true_corum_gecko.pdf"),
           width=6, height=4)

corum_wang_true <-
    run_int_ext_true(corum_wang_dat$edgelist,
                     corum_wang_dat$genesets,
                     rank_thresholds)

corum_wang_true %>%
    edgedens_true_gg(title="CORUM mapped onto Wang2017 Similarity Network") +
    ggsave(file.path(out_dir, "edgedens_true_corum_wang.pdf"),
           width=6, height=4)


corum_true <- bind_rows(
    corum_avana_true %>% mutate(Network="Avana"),
    corum_avana_2017_true %>% mutate(Network="Avana 2017"),
    corum_rnai_true %>% mutate(Network="RNAi"),
    corum_coxpr_true %>% mutate(Network="COXPRESdb"),
    corum_gecko_true %>% mutate(Network="GeCKOv2"),
    corum_wang_true %>% mutate(Network="Wang 2017")
)


corum_shuffle <- bind_rows(
    corum_avana_shuffle %>% mutate(Network="Avana"),
    corum_avana_2017_shuffle %>% mutate(Network="Avana 2017"),
    corum_rnai_shuffle %>% mutate(Network="RNAi"),
    corum_coxpr_shuffle %>% mutate(Network="COXPRESdb"),
    corum_gecko_shuffle %>% mutate(Network="GeCKOv2"),
    corum_wang_shuffle %>% mutate(Network="Wang 2017")
)



empirical_shuffle_results <-
    inner_join(corum_true, corum_shuffle) %>%
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
    mutate(Recall = cumsum(BestFDR < 0.05)/length(corum_list))

write_tsv(rank_fdrs, file.path(out_dir, "corum_shuffle_10k_results.tsv"))


ggplot(rank_fdrs, aes(FirstRank, Recall, color=Network)) +
    geom_hline(aes(yintercept=RecallFinal, color=Network), linetype=2,
               data=rank_fdrs %>% group_by(Network) %>%
                   summarize(RecallFinal = sum(BestFDR<0.05) /
                                 length(corum_list))) +
    geom_line() +
    scale_y_continuous(limits=c(0, 0.5), expand=c(0, 0)) +
    scale_x_continuous(trans="log2", breaks=rank_thresholds) +
    labs(x="Rank Threshold", y="Recall of Complexes",
         title="CORUM mapped onto Similarity Networks") +
    ggsave(file.path(out_dir, "rank_recall_corum.pdf"),
           width=6, height=4)

rank_fdrs %>%
    group_by(Network) %>%
    arrange(BestFDR) %>%
    mutate(Recall = cumsum(BestFDR < 1)/length(corum_list)) %>%
    ggplot(aes(Recall, 1-BestFDR, color=Network)) +
    geom_hline(yintercept=0.95, linetype=2) +
    geom_vline(aes(xintercept=RecallFinal, color=Network), linetype=2,
               data=rank_fdrs %>% group_by(Network) %>%
                   summarize(RecallFinal = sum(BestFDR<0.05) /
                                 length(corum_list))) +
    geom_line() +
    scale_y_continuous(limits=c(0.5, 1.02), expand=c(0, 0)) +
    scale_x_continuous(expand=c(0, 0)) +
    labs(y="Precision", x="Recall",
         title="CORUM Precision-Recall in Similarity Networks") +
    ggsave(file.path(out_dir, "precision_recall_corum.pdf"),
           width=6, height=4)


empirical_shuffle_results %>%
    filter(Network != "Gecko") %>%
    group_by(Network, Geneset) %>%
    arrange(Rank) %>%
    mutate(Right = c(FDR[2:n()], NA),
           Left = c(NA, FDR[1:(n()-1)])) %>%
    filter(FDR < 0.05) %>%
    select(Network, Geneset, Rank, FDR, Left, Right) %>%
    gather(Group, Sig, Left, FDR, Right) %>%
    mutate(Group = factor(Group, levels=c("Left", "FDR", "Right"), ordered=T)) %>%
    ggplot(aes(x=as.factor(Rank), y=-log10(Sig), color=Group)) +
    geom_hline(yintercept=-log10(0.05)) +
    geom_boxplot( position=position_dodge(width=0.9)) +
    geom_text(aes(x=as.factor(Rank), label=NSig),
              y=max(-log10(empirical_shuffle_results$FDR), na.rm=T),
              vjust=-0.05,
              data=empirical_shuffle_results %>%
                  filter(Network != "Gecko") %>%
                  group_by(Network, Rank) %>%
                  summarise(NSig = sum(FDR < 0.05, na.rm=T)),
              inherit.aes=F) +
    facet_wrap(~ Network, ncol=2, scales="free") +
    scale_color_manual(values=brewer_pal(palette="Set2")(2) %>%
                           set_names(c("Left", "Right")) %>%
                           c(FDR="Black")) +
    theme(strip.background=element_blank()) +
    ggsave(file.path(out_dir, "sig_levels_across_ranks.pdf"),
           width=10, height=10)



near_misses <-
    empirical_shuffle_results %>%
    group_by(Network, Geneset) %>%
    arrange(Rank) %>%
    summarise(BestSig = min(FDR, na.rm=T),
              RankOfBestSig = Rank[FDR == BestSig][1],
              NearestMiss = min(FDR[FDR > 0.05], na.rm=T),
              RankOfNearestMiss = Rank[FDR == NearestMiss][1])

near_misses %>% filter(!is.infinite(BestSig),
                       !is.infinite(NearestMiss)) %>%
    ggplot(aes(-log10(NearestMiss), -log10(BestSig),
               color=1/RankOfNearestMiss)) +
    facet_wrap(~ Network, scales="free", ncol=3) +
    geom_hline(yintercept=-log10(0.05)) + geom_vline(xintercept=-log10(0.05)) +
    geom_point(position=position_jitter(width=0.05,
                                        height=0.05),
               show.legend=F) +
    scale_color_continuous(trans="log2") +
    ggsave(file.path(out_dir, "near_fdr_misses_corum.pdf"),
           width=5, height=5)

