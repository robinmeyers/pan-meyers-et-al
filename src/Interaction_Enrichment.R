library(ProjectTemplate)
load.project()

if (config$threads > 1) {
    library(doMC)
    registerDoMC(cores=config$threads)
    do_parallel <- T
} else {
    do_parallel <- F
}

plot_dir <- file.path("./graphs/interaction_enrichment", Sys.Date())
dir.create(plot_dir, recursive = T, showWarnings = F)

corum_pairs <- data_frame(Complex = names(corum_list),
                          Gene = corum_list) %>%
    unnest(Gene) %>%
    inner_join(., ., by="Complex") %>%
    filter(Gene.x != Gene.y) %>%
    alphabetize_genes()

humap_pairs <- data_frame(Complex = seq_along(humap_list),
                          Gene = humap_list) %>%
    unnest(Gene) %>%
    inner_join(., ., by="Complex") %>%
    filter(Gene.x != Gene.y) %>%
    alphabetize_genes()

genepair_lists <- list(CORUM = corum_pairs,
                       hu.MAP = humap_pairs,
                       Y2H = rolland2014_interactome_y2h,
                       Subcell = thul2017_subcellular_y2h,
                       Mann = hein2015_interactome_qubic,
                       BioPlex2 = huttlin2017_bioplex_ppi,
                       Marcottte = wan2015_complexes_cofrac)

rnai_corr_interaction_enrichment <-
    ldply(genepair_lists, function(gp) {
        interaction_enrichment(rnai_dep_corr, genepairs = gp,
                                    binsize=1000, numbins=100)
    }, .id="Source", .parallel = do_parallel)

avana_corr_interaction_enrichment <-
    ldply(genepair_lists, function(gp) {
        interaction_enrichment(avana_dep_corr, genepairs = gp,
                               binsize=1000, numbins=100)
    }, .id="Source", .parallel = do_parallel)


rnai_gg <- ggplot(rnai_corr_interaction_enrichment,
                 aes(x=CoDep, y=log2(estimate), color=Source)) +
    geom_hline(yintercept=0, size=0) +
    geom_line(size=1.5) +
    scale_color_brewer(palette = "Paired") +
    scale_y_continuous(expand = c(0, 0),
                       breaks = pretty_breaks(4)) +
    labs(x="Co-dependency Overall Rank",
         y="Fold Enrichment (log2)",
         title="RNAi") +
    theme(legend.position="right",
          legend.title = element_blank())

crispr_gg <- ggplot(avana_corr_interaction_enrichment,
            aes(x=CoDep, y=log2(estimate), color=Source)) +
    geom_hline(yintercept=0, size=0) +
    geom_line(size=1.5) +
    scale_color_brewer(palette = "Paired") +
    scale_y_continuous(expand = c(0, 0),
                       breaks = pretty_breaks(4)) +
    labs(x="Co-dependency Overall Rank",
         y="Fold Enrichment (log2)",
         title="CRISPR") +
    theme(legend.position="right",
          legend.title = element_blank())


rnai_gg +
    ggsave(file.path(plot_dir, "rnai_interaction_enrichment.pdf"),
           width=6, height=4)

crispr_gg +
    ggsave(file.path(plot_dir, "crispr_interaction_enrichment.pdf"),
           width=6, height=4)

