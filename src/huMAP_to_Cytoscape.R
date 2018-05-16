library(ProjectTemplate)
load.project(override.config = list(cache_loading=F, munging=F))

load("./cache/humap_list.RData")
load("./cache/corum_list.RData")
load("./cache/avana_2017_dep_rank.RData")


#humap complexes export to Cytoscape
humap_thresholds <- read_tsv("./data/interim/humap_thresholds.tsv") %>%
    mutate(Geneset = Geneset %>% as.character)
humap_shuffle_results <- read_tsv("./data/interim/humap_shuffle_10k_results.tsv") %>%
    mutate(Geneset = Geneset %>% as.character)

humap_pairs <- humap_list %>% list_to_pairs(unique_pairs=F)

crispr_ranked_df <- avana_2017_dep_rank %>%
    mat.to.df("Gene.x", "Gene.y", "Rank") %>%
    filter(Gene.x < Gene.y)

humap_crispr <- humap_thresholds %>%
    select(Geneset, Threshold = Rank) %>%
    left_join(humap_pairs, by = c("Geneset" = "Complex")) %>%
    inner_join(crispr_ranked_df) %>%
    filter(Rank <= Threshold)

write_tsv(humap_crispr,
          "./data/interim/cytoscape_humap_crispr_nets.txt")

#write humap/corum gene annotation
humap_genes <- data_frame(Gene = unlist(humap_list) %>% unique) %>%
    mutate(In_CORUM = Gene %in% unlist(corum_list))

write_tsv(humap_genes,
          "./data/interim/humap_genes_in_corum.txt")

