library(ProjectTemplate)
load.project(override.config = list(cache_loading=F, munging=F))

### assign complexes to RNAi or CRISPR for downstream analysis
load("./cache/corum_list.RData")
load("./cache/avana_2017_dep_rank.RData")
load("./cache/rnai_dep_rank.RData")

out_dir <- file.path("./output/complex_assign/", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)

if (config$threads > 1) {
    library(doMC)
    registerDoMC(cores=config$threads)
    do_parallel <- T
} else {
    do_parallel <- F
}

#CORUM results and edgeweights
corum_thresholds <- read_tsv("./data/interim/corum_thresholds.tsv") %>%
    filter(Network %in% c("Avana 2017", "RNAi"))
corum_shuffle_results <- read_tsv("./data/interim/corum_shuffle_10k_results.tsv") %>%
    filter(Network %in% c("Avana 2017", "RNAi"))

signif_complexes <- corum_shuffle_results %>%
    filter(BestFDR < 0.05) %>%
    pull(Geneset) %>% unique

rnai_max_density <- rnai_dep_rank %>%
    calculate_graph_thresholds(corum_list[signif_complexes],
                               metric="ratio", do_parallel=T)
crispr_max_density <- avana_2017_dep_rank %>%
    calculate_graph_thresholds(corum_list[signif_complexes],
                               metric="ratio", do_parallel=T)

#Load CRISPR vs RNAi assignments
shuffle_assignment <- corum_shuffle_results %>%
    filter(BestFDR < 0.05) %>%
    group_by(Geneset) %>%
    summarise(Dataset = case_when(
        any(Network=="Avana 2017") & any(Network=="RNAi") ~ "Both",
        any(Network=="Avana 2017") ~ "CRISPR",
        any(Network=="RNAi") ~ "RNAi")) %>%
    mutate(Dataset = factor(Dataset, levels = c("RNAi", "Both", "CRISPR")))



max_densities <- full_join(rnai_max_density,
                        crispr_max_density,
                        by = "Geneset") %>%
    inner_join(shuffle_assignment)

#complex sorting plot
max_densities %>%
    select(Geneset, Max_RNAi = Max.x, Max_CRISPR = Max.y, Dataset) %>%
    ggplot(aes(Max_CRISPR,
               Max_RNAi,
               color = Dataset,
               label = Geneset)) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    geom_abline(slope=1) +
    scale_color_manual(values = c("#00BFC4", "#3B4CA0", "#F8766D")) +
    labs(x = "Maximum Edge Density Ratio (CRISPR)",
         y = "Maximum Edge Density Ratio (RNAi)") +
    theme(legend.position = "bottom") +
    ggsave(file.path(out_dir, "edge_density_scatterplot.pdf"),
           width = 4.25, height = 4.5)


crispr_ranked_df <- avana_2017_dep_rank %>%
    mat.to.df("Gene.x", "Gene.y", "Rank") %>%
    filter(Gene.x < Gene.y,
           Gene.x %in% unlist(corum_list),
           Gene.y %in% unlist(corum_list))

rnai_ranked_df <- rnai_dep_rank %>%
    mat.to.df("Gene.x", "Gene.y", "Rank") %>%
    filter(Gene.x < Gene.y,
           Gene.x %in% unlist(corum_list),
           Gene.y %in% unlist(corum_list))

corum_pairs <- corum_list %>%
    list_to_pairs(unique_pairs = F)

#sort into either CRISPR or RNAi, not Both
dataset_assignment <- max_densities %>%
    select(Geneset, Max_RNAi = Max.x, Max_CRISPR = Max.y, Dataset) %>%
    filter(!is.na(Dataset)) %>%
    mutate(Binary = case_when(is.na(Max_RNAi) ~ "CRISPR",
                              is.na(Max_CRISPR) ~ "RNAi",
                              Max_RNAi > Max_CRISPR ~ "RNAi",
                              Max_CRISPR > Max_RNAi ~ "CRISPR"))


corum_crispr_edges <- dataset_assignment %>%
    filter(Binary == "CRISPR") %>%
    select(Geneset) %>%
    left_join(corum_thresholds %>%
                  filter(Network == "Avana 2017") %>%
                  select(Geneset, Threshold = Rank)) %>%
    left_join(corum_pairs, by=c("Geneset"="Complex")) %>%
    inner_join(crispr_ranked_df) %>%
    filter(Rank <= Threshold)

corum_rnai_edges <- dataset_assignment %>%
    filter(Binary == "RNAi") %>%
    select(Geneset) %>%
    left_join(corum_thresholds %>%
                  filter(Network == "RNAi") %>%
                  select(Geneset, Threshold = Rank)) %>%
    left_join(corum_pairs, by=c("Geneset"="Complex")) %>%
    inner_join(rnai_ranked_df) %>%
    filter(Rank <= Threshold)

corum_crispr_edges %>% select(Gene.x, Gene.y) %>% distinct() %>%
    write_tsv(file.path(out_dir, "cytoscape_corum_crispr_nets.tsv"))
corum_rnai_edges %>% select(Gene.x, Gene.y) %>% distinct() %>%
    write_tsv(file.path(out_dir, "cytoscape_corum_rnai_nets.tsv"))
