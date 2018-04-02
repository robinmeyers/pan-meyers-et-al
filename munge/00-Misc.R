
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


cache("ccds", {
    read_tsv("./data/raw/CCDS.current.txt",
             col_types = cols(.default = col_character())) %>%
        mutate(cds_from = as.integer(cds_from),
               cds_to = as.integer(cds_to))
})
