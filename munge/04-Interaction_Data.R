
cache("huttlin2017_bioplex_ppi", {
    read_tsv("./data/raw/nature22366-s2.tsv",
             col_types = "----cc--n") %>%
        set_colnames(c("Gene.x", "Gene.y", "Confidence")) %>%
        alphabetize_genes(Gene.x, Gene.y)
})

cache("thul2017_subcellular_y2h", {
    read_tsv("./data/raw/aal3321_Thul_SM_table_S17.tsv",
             col_types = "cc--c") %>%
        set_colnames(c("Gene.x", "Gene.y", "Localization")) %>%
        de_excelify_genes(Gene.x) %>%
        de_excelify_genes(Gene.y) %>%
        alphabetize_genes(Gene.x, Gene.y)
})

cache("rolland2014_interactome_y2h", {
    read_tsv("./data/raw/rolland2014_mmc3.tsv",
             col_types = "--cc-") %>%
        set_colnames(c("Gene.x", "Gene.y")) %>%
        de_excelify_genes(Gene.x) %>%
        de_excelify_genes(Gene.y) %>%
        alphabetize_genes(Gene.x, Gene.y)
})


cache("wan2015_complexes_cofrac", {
    raw <- read_tsv("./data/raw/wan2015_table_S4.tsv",
                         col_types = 'cicc') %>%
        mutate(Gene = str_split(GeneName, ";")) %>%
        select(ComplexID, Gene) %>%
        unnest()

    inner_join(raw, raw, by="ComplexID") %>%
        filter(Gene.x < Gene.y) %>%
        alphabetize_genes()
})


cache("hein2015_interactome_qubic", {
    mann <- read_tsv("./data/raw/hein2015_mmc3.tsv") %>%
        select(Gene.x = bait.Gene.name, Gene.y = prey.Gene.name) %>%
        alphabetize_genes()
})
