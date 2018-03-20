
cache("ccle_rpkm", {
    read_csv("./data/raw/CCLE_RNAseq_081117.rpkm.csv",
             col_types = cols(X1 = col_character(),
                              .default=col_double())) %>%
        as.data.frame() %>%
        set_rownames(.[,1]) %>%
        select(-1) %>%
        as.matrix() %>%
        log2
})

cache("hart_nonessentials", {
    hart_ctrls <- read_tsv("./data/raw/hart2014_seed_core_essentials.txt")
    ne_genes <- hart_ctrls$`Nonessential Genes (NE)` %>% str_subset("[A-Z0-9]+")
})

cache("hart_coreessentials", {
    hart_ctrls <- read_tsv("./data/raw/hart2014_seed_core_essentials.txt")
    cce_genes <- hart_ctrls$`ConstitutiveCoreEssentials(CCE)` %>% str_subset("[A-Z0-9]+")
})
