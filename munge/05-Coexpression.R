cache("coxpres_db", {
    raw <- read_csv("./data/raw/coexpressdb.csv",
                    col_types = cols(X1 = col_character(),
                                     .default=col_double()))

    entrez_to_hugo <- ccds %>%
        distinct(gene, gene_id) %$%
        set_names(gene, gene_id)

    raw %>% as.data.frame %>%
        set_rownames(.$X1) %>%
        select(-X1) %>%
        as.matrix() %>%
        set_rownames(entrez_to_hugo[rownames(.)]) %>%
        set_colnames(entrez_to_hugo[colnames(.)])
}, depends = "ccds")
