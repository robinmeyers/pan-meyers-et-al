edge_threshold <- function(geneset, edge_df, metric="difference") {
    ranked_df <- edge_df %>%
        filter(Gene.x %in% geneset | Gene.y %in% geneset) %>%
        mutate(In_Complex = Gene.x %in% geneset & Gene.y %in% geneset) %>%
        group_by(Edgeweight) %>%
        summarize(Internal = sum(In_Complex),
                  External = sum(!In_Complex)) %>%
        ungroup %>%
        arrange(desc(Edgeweight)) %>%
        mutate(Int_Edge_Density = cumsum(Internal)/sum(Internal),
               Ext_Edge_Density = cumsum(External)/sum(External),
               Int_Edge_Density_Weighted = cumsum(Internal*Edgeweight)/sum(Internal*Edgeweight),
               Ext_Edge_Density_Weighted = cumsum(External*Edgeweight)/sum(External*Edgeweight),
               Diff_Edge_Density = Int_Edge_Density_Weighted - Ext_Edge_Density_Weighted,
               Edge_Density_Ratio = Int_Edge_Density / Ext_Edge_Density)

    if (metric == "difference") {
        output <- ranked_df %>%
            arrange(desc(Diff_Edge_Density), Edgeweight) %>%
            slice(1) %>%
            mutate(Max = Diff_Edge_Density)
    } else if (metric == "ratio") {
        output <- ranked_df %>%
            arrange(desc(Edge_Density_Ratio), Edgeweight) %>%
            slice(1) %>%
            mutate(Max = Edge_Density_Ratio)
    } else {
        stop("Invalid metric argument")
    }

    return(output)
}


calculate_graph_thresholds <- function(mat, geneset_list,
                                       metric="difference", do_parallel = FALSE) {

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
        ldply(edge_threshold, edge_df=rank_df, metric=metric,
              .parallel=do_parallel, .id="Geneset") %>%
        transmute(Geneset = as.character(Geneset),
                  Rank = 1/Edgeweight,
                  Max = Max)
}

