edge_threshold <- function(geneset, edge_df) {
    ranked_df <- edge_df %>%
        filter(Gene.x %in% geneset | Gene.y %in% geneset) %>%
        mutate(In_Complex = Gene.x %in% geneset & Gene.y %in% geneset) %>%
        group_by(Edgeweight) %>%
        summarize(Internal = sum(In_Complex),
                  External = sum(!In_Complex)) %>%
        ungroup %>%
        arrange(desc(Edgeweight)) %>%
        mutate(Int_Edge_Density = cumsum(Internal*Edgeweight)/sum(Internal*Edgeweight),
               Ext_Edge_Density = cumsum(External*Edgeweight)/sum(External*Edgeweight),
               Diff_Edge_Density = Int_Edge_Density - Ext_Edge_Density)

    output <- ranked_df %>% arrange(desc(Diff_Edge_Density), Edgeweight) %>% slice(1)
    return(output)
}


calculate_graph_thresholds <- function(mat, geneset_list, do_parallel = FALSE) {

    rank_df <-
        mat %>%
        mat.to.df("Gene.x", "Gene.y", "Rank") %>%
        filter(Gene.x < Gene.y) %>%
        mutate(Edgeweight = 1/Rank)

    # remove complexes that have fewer than 2 members in the screening dataset
    genesets_to_test <- llply(geneset_list, function(genes) {
        screened <- genes[genes %in% rownames(mat)] %>% unique()
        if (length(screened) < 2)
            return(NA)
        else
            return(screened)
    }) %>% magrittr::extract(!is.na(.))

    genesets_to_test %>%
        ldply(edge_threshold, edge_df=rank_df,
              .parallel=do_parallel, .id="Geneset") %>%
        transmute(Geneset = as.character(Geneset),
                  Rank = 1/Edgeweight)
}
