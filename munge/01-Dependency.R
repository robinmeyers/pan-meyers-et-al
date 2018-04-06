cache("avana_dep", {
    read_csv("./data/raw/avana_18Q1_gene_effects.csv",
             col_types = cols(X1 = col_character(),
                              .default=col_double())) %>%
        as.data.frame() %>%
        set_rownames(.[,1]) %>%
        select(-1) %>%
        as.matrix()
})

cache("avana_2017_dep", {
    read_csv("./data/raw/avana_2017_gene_effects.csv",
             col_types = cols(X1 = col_character(),
                              .default=col_double())) %>%
        as.data.frame() %>%
        set_rownames(.[,1]) %>%
        select(-1) %>%
        as.matrix()
})

cache("rnai_dep", {
    read_csv("./data/raw/rnai_gene_effects.csv",
             col_types = cols(X1 = col_character(),
                              .default=col_double())) %>%
        as.data.frame() %>%
        set_rownames(.[,1]) %>%
        select(-1) %>%
        as.matrix()
})

cache("gecko_dep", {
    read_csv("./data/raw/gecko_gene_effects.csv",
             col_types = cols(X1 = col_character(),
                              .default=col_double())) %>%
        as.data.frame() %>%
        set_rownames(.[,1]) %>%
        select(-1) %>%
        as.matrix()
})

cache("wang_dep", {
    read_csv("./data/raw/wang_gene_effects.csv",
             col_types = cols(X1 = col_character(),
                              .default=col_double())) %>%
        as.data.frame() %>%
        set_rownames(.[,1]) %>%
        select(-1) %>%
        as.matrix()
})


cache("avana_analysis_genes", {
    globalSD <- sd(avana_dep, na.rm=T)
    most_dependent_line <- adply(avana_dep, 1, function(x) {
        xmin <- which.min(x)
        data_frame(CellLine = names(xmin),
                   Dep = x[xmin],
                   SD = ((x[xmin] - mean(x, na.rm=T)) / globalSD))
    }, .id="Gene")

    xpr <- most_dependent_line %>%
        select(CellLine, Gene) %>%
        filter(CellLine %in% colnames(ccle_rpkm),
               Gene %in% rownames(ccle_rpkm)) %>%
        mutate(RPKM = ccle_rpkm[as.matrix(.[,c("Gene", "CellLine")])])

    analysis_genes <- most_dependent_line %>%
        left_join(xpr, by=c("Gene", "CellLine")) %>%
        filter(Dep < -0.3 & (is.na(RPKM) | RPKM > 0))
    return(as.character(analysis_genes$Gene))
}, depends = c("avana_dep", "ccle_rpkm"))



cache("avana_2017_analysis_genes", {
    globalSD <- sd(avana_2017_dep, na.rm=T)
    most_dependent_line <- adply(avana_2017_dep, 1, function(x) {
        xmin <- which.min(x)
        data_frame(CellLine = names(xmin),
                   Dep = x[xmin],
                   SD = ((x[xmin] - mean(x, na.rm=T)) / globalSD))
    }, .id="Gene")

    xpr <- most_dependent_line %>%
        select(CellLine, Gene) %>%
        filter(CellLine %in% colnames(ccle_rpkm),
               Gene %in% rownames(ccle_rpkm)) %>%
        mutate(RPKM = ccle_rpkm[as.matrix(.[,c("Gene", "CellLine")])])

    analysis_genes <- most_dependent_line %>%
        left_join(xpr, by=c("Gene", "CellLine")) %>%
        filter(Dep < -0.3 & (is.na(RPKM) | RPKM > 0))
    return(as.character(analysis_genes$Gene))
}, depends = c("avana_2017_dep", "ccle_rpkm"))

cache("rnai_analysis_genes", {
    read_tsv("./data/raw/rnai_analysis_genes.tsv",
             col_types = cols(gene = col_character()))$gene
})
