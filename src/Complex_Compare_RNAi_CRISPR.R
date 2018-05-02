
library(ggsignif)

library(ProjectTemplate)
load.project(override.config = list(cache_loading=F, munging=F))

out_dir <- file.path("./output/complex_compare", Sys.Date())
dir.create(out_dir, recursive = T, showWarnings = F)

load("./cache/corum_list.RData")
load("./cache/ccle_rpkm.RData")

corum_significant_genesets <- read_tsv("./data/interim/corum_shuffle_10k_sig_genesets.tsv")

dnds <- read_tsv("./data/raw/human_mouse_dnds.txt") %>%
    filter(!is.na(`dN with Mouse`), `dS with Mouse` > 0) %>%
    mutate(dNdS = `dN with Mouse`/`dS with Mouse`) %>%
    select(Gene = `Gene name`, dNdS) %>%
    distinct()

human_tax_id <- 9606
phy_retent <- read_tsv("./data/raw/homologene.data.txt",
                       col_names = c("HID", "Taxonomy_ID",
                                     "Gene_ID", "Gene_Symbol",
                                     "Protein_gi", "Protein_accession")) %>%
    group_by(HID) %>%
    summarize(Gene = Gene_Symbol[Taxonomy_ID == human_tax_id][1],
              Num_Phyla = length(unique(Taxonomy_ID))) %>%
    filter(!is.na(Gene))

ccle_rpkm_75 <- adply(ccle_rpkm, 1, quantile, probs=0.75, na.rm = TRUE) %>%
    set_colnames(c("Gene", "RPKM"))

ceg1 <- read_tsv("./data/raw/hart2014_seed_core_essentials.txt")$`ConstitutiveCoreEssentials(CCE)` %>%
    magrittr::extract(!is.na(.))
ceg2 <- read_tsv("./data/raw/CEGv2.txt")$Gene
core_ess <- union(ceg1, ceg2)

corum_stats <- data_frame(Geneset = names(corum_list),
                          Gene = corum_list) %>%
    unnest(Gene) %>%
    left_join(dnds) %>%
    left_join(phy_retent) %>%
    left_join(ccle_rpkm_75) %>%
    group_by(Geneset) %>%
    summarize(Median_dNdS = median(dNdS, na.rm = TRUE),
              Median_Expr = median(RPKM, na.rm = TRUE),
              Frac_Essential = sum(Gene %in% core_ess)/n())

significant_complexes <- read_tsv("./data/interim/corum_shuffle_10k_sig_genesets.tsv") %>%
    filter(Network %in% c("Avana 2017", "RNAi")) %>%
    group_by(Geneset) %>%
    summarise(Dataset = case_when(
        any(Network=="Avana 2017") & any(Network=="RNAi") ~ "Both",
        any(Network=="Avana 2017") ~ "CRISPR",
        any(Network=="RNAi") ~ "RNAi")) %>%
    mutate(Dataset = factor(Dataset, levels = c("RNAi", "Both", "CRISPR")))


stat_name <- c("Frac_Essential","Median_Expr", "Median_dNdS")
stat_labels <- c("Fraction of Essential Genes","Median Gene Expression", "Median dNdS")
complex_compare <- left_join(significant_complexes, corum_stats,
                              by = "Geneset") %>%
    gather(Stat, Value, -Geneset, -Dataset) %>%
    filter(Stat %in% stat_name) %>%
    mutate(Stat = factor(Stat, levels = stat_name, labels = stat_labels))

ggplot(complex_compare, aes(Dataset, Value, fill=Dataset)) +
    geom_boxplot() +
    facet_wrap(~ Stat, scales = "free", dir = "v", nrow = 1) +
    geom_signif(comparisons = list(c("RNAi", "Both"),
                                   c("Both", "CRISPR"),
                                   c("CRISPR", "RNAi")),
                step_increase = c(0, 0.05, 0.05), map_signif_level = TRUE) +
    theme_minimal() +
    scale_fill_manual(values = c("#00BFC4", "#3B4CA0","#F8766D")) +
    ggsave(file.path(out_dir, "complex_compare_boxplot.pdf"),
       width = 8, height = 4)

