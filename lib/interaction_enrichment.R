interaction_enrichment <- function(mat, genepairs,
                                   binsize=1000, numbins=100) {
    codep_df <- mat %>%
        magrittr::extract(rownames(.) %in% c(genepairs$Gene.x, genepairs$Gene.y),
                          colnames(.) %in% c(genepairs$Gene.x, genepairs$Gene.y)) %>%
        mat.to.df("Gene.x", "Gene.y", "Edgeweight") %>%
        filter(Gene.x < Gene.y)


    genepairs <- genepairs %>%
        distinct(Gene.x, Gene.y) %>%
        filter(Gene.x %in% rownames(mat),
               Gene.y %in% rownames(mat)) %>%
        mutate(Interaction = T)

    total_pairs <- nrow(codep_df)
    total_interacting_pairs <- nrow(genepairs)

    combined_df <- left_join(codep_df, genepairs, by=c("Gene.x", "Gene.y")) %>%
        mutate(Interaction = !is.na(Interaction),
               EW_Rank = rank(-Edgeweight, ties.method="min"))


    ranked_df <- combined_df %>%
        mutate(EW_Bin = cut(EW_Rank, breaks=seq(0, total_pairs+binsize, by=binsize),
                            ordered_result=T)) %>%
        group_by(EW_Bin) %>%
        summarise(CoDepBin = n(),
                  InteractBin = sum(Interaction)) %>%
        ungroup() %>%
        arrange(EW_Bin) %>%
        head(n=numbins) %>%
        mutate(CoDep = cumsum(CoDepBin),
               CoDep_Interact = cumsum(InteractBin),
               CoDep_NoInteract = CoDep - CoDep_Interact,
               NoCoDep_Interact = total_interacting_pairs - CoDep_Interact,
               NoCoDep_NoInteract = total_pairs - CoDep - NoCoDep_Interact)


    fisher_df <- ranked_df %>%
        rowwise() %>%
        do(broom::glance(fisher.test(
            matrix(c(.$CoDep_Interact, .$CoDep_NoInteract,
                     .$NoCoDep_Interact, .$NoCoDep_NoInteract),
                   nrow=2))))


    return(cbind(ranked_df, fisher_df))
}
