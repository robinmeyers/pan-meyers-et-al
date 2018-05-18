library(ProjectTemplate)
load.project(override.config = list(cache_loading=F, munging=F))

if (config$threads > 1) {
    library(doMC)
    registerDoMC(cores=config$threads)
    do_parallel <- T
} else {
    do_parallel <- F
}

load("./cache/corum_list.RData")
load("./cache/humap_list.RData")

load("./cache/avana_dep_rank.RData")
load("./cache/avana_2017_dep_rank.RData")
load("./cache/rnai_dep_rank.RData")
# load("./cache/gecko_dep_rank.RData")
# load("./cache/wang_dep_rank.RData")
load("./cache/coxpres_db.RData")
coxpres_db_rank <- edgeweight_symmetric_rank(coxpres_db)

corum_significant_genesets <- read_tsv("./data/interim/corum_shuffle_10k_sig_genesets.tsv")
humap_significant_genesets <- read_tsv("./data/interim/humap_shuffle_10k_sig_genesets.tsv")


avana_corum_thresholds <-
    calculate_graph_thresholds(avana_dep_rank,
                               corum_list[corum_significant_genesets %>%
                                              filter(Network == "Avana") %>%
                                              pull(Geneset)],
                               do_parallel = T)

avana_2017_corum_thresholds <-
    calculate_graph_thresholds(avana_2017_dep_rank,
                               corum_list[corum_significant_genesets %>%
                                              filter(Network == "Avana 2017") %>%
                                              pull(Geneset)],
                               do_parallel = T)

rnai_corum_thresholds <-
    calculate_graph_thresholds(rnai_dep_rank,
                               corum_list[corum_significant_genesets %>%
                                              filter(Network == "RNAi") %>%
                                              pull(Geneset)],
                               do_parallel = T)

coxpresdb_corum_thresholds <-
    calculate_graph_thresholds(coxpres_db_rank,
                               corum_list[corum_significant_genesets %>%
                                              filter(Network == "COXPRESdb") %>%
                                              pull(Geneset)],
                               do_parallel = T)

corum_thresholds <-
    bind_rows(avana_corum_thresholds %>% mutate(Network = "Avana"),
              avana_2017_corum_thresholds %>% mutate(Network = "Avana 2017"),
              rnai_corum_thresholds %>% mutate(Network = "RNAi"),
              coxpresdb_corum_thresholds %>% mutate(Network = "COXPRESdb"))


corum_thresholds %>%
    select(Network, Geneset, Rank) %>%
    write_tsv("./data/interim/corum_thresholds.tsv")


avana_2017_humap_thresholds <-
    calculate_graph_thresholds(avana_2017_dep_rank,
                               humap_list[humap_significant_genesets %>%
                                              filter(Network == "Avana 2017") %>%
                                              pull(Geneset)],
                               do_parallel = T)

humap_thresholds <-
    bind_rows(#avana_humap_thresholds %>% mutate(Network = "Avana"),
              avana_2017_humap_thresholds %>% mutate(Network = "Avana 2017")
              # rnai_humap_thresholds %>% mutate(Network = "RNAi")
              )


humap_thresholds %>%
    select(Network, Geneset, Rank) %>%
    write_tsv("./data/interim/humap_thresholds.tsv")


