library(ProjectTemplate)
load.project(override.config = list(cache_loading=F, munging=F))

#Correlation distributions for CORUM complex
load("./cache/corum_list.RData")
load("./cache/avana_2017_dep_corr.RData")
load("./cache/rnai_dep_corr.RData")

out_dir <- file.path("./output/correlation_distributions", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=T)

corum_pairs <- list_to_pairs(corum_list, unique_pairs=F)
corum_genes <- corum_list %>% unlist() %>% unique()

#Filter corum genes by presence in fitness dataset
labeled_corr_df <- function(mat, genes) {
    genes <- intersect(genes, rownames(mat))
    corr_df <- mat.to.df(mat[genes, genes],
                         "Gene.x", "Gene.y", "Corr") %>%
        filter(Gene.x < Gene.y) %>%
        left_join(corum_pairs) %>%
        mutate(Same_Complex = !is.na(Complex)) %>%
        distinct(Gene.x, Gene.y, Corr, Same_Complex)
}

corr_density <- function(corr_df) {
    corr_df %>%
        ggplot(aes(x = Corr, color = Same_Complex)) +
        geom_density() +
        xlim(-1, 1) +
        labs(x = "Pearson Correlation",
             y = "Density") +
        scale_color_manual(values = c("#387EB9", "#E31F26")) +
        theme(legend.position="none")
}

corum_corr_dist <- function(cor) {
    labeled_corr_df(cor, corum_genes) %>% corr_density()
}

#Adding statistics

calculate_ks_pval_for_corum <- function(corr) {
    df  <- labeled_corr_df(corr, corum_genes)

    ks.test(df %>% filter(Same_Complex == F) %>% pull(Corr),
            df %>% filter(Same_Complex == T) %>% pull(Corr))
}

calculate_ks_pval_for_corum(avana_2017_dep_corr)
calculate_ks_pval_for_corum(rnai_dep_corr)

#Correlation distributions across CORUM
plot_grid(plotlist = list(corum_corr_dist(rnai_dep_corr) + ggtitle("RNAi"),
                          corum_corr_dist(avana_2017_dep_corr) + ggtitle("CRISPR"))) +
    ggsave(file.path(out_dir, "correlation_distributions.pdf"),
       width = 10, height = 5)
