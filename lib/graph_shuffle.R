

prepare_dataset <- function(corr_mat, full_genesets) {

    genesets <- llply(full_genesets, function(s) {
        intersect(s, rownames(corr_mat))
    })

    geneset_df <- data_frame(Geneset=names(genesets),
                             Gene = genesets) %>%
        unnest(Gene)


    rank_mat <- corr_mat %>%
        edgeweight_symmetric_rank() %>%
        magrittr::extract(rownames(.) %in% geneset_df$Gene,
                          colnames(.) %in% geneset_df$Gene)

    edgelist <-
        rank_mat %>%
        mat.to.df("Gene.x", "Gene.y", "Rank") %>%
        filter(Gene.x < Gene.y)

    return(list(edgelist=edgelist, genesets=geneset_df))
}


run_int_ext_true <- function(edgelist, geneset, rank_thresholds) {
    ldply(rank_thresholds, function(rank_threshold) {

        cat("Running on true graph at", rank_threshold, "\n")

        g <- edgelist %>% filter(Rank <= rank_threshold) %>%
            select(from=Gene.x, to=Gene.y) %>%
            graph_from_data_frame(directed=F,
                                  vertices=data_frame(Gene = unique(geneset$Gene)))

        int_ext_ratio <- calc_int_ext_dens_ratio(g, geneset, metric="ratio")

        data_frame(Rank=rank_threshold,
                   Geneset = names(int_ext_ratio),
                   IntExtTrue = int_ext_ratio)
    })
}


edgedens_true_gg <- function(int_ext_true,
                             title=NULL) {
    int_ext_true %>%
        mutate(IntExtTrue = ifelse(IntExtTrue==0,
                                   min(IntExtTrue[IntExtTrue>0], na.rm=T)/10,
                                   IntExtTrue)) %>%
        ggplot(aes(x=IntExtTrue, y=as.factor(Rank), fill=1/Rank)) +
        geom_density_ridges(scale=5, show.legend=F) +
        scale_x_continuous(expand=c(0,0), trans="log10", breaks=c(1, 10, 100, 1000)) +
        scale_fill_continuous(trans="log2") +
        labs(x="Internal-to-External Edge Density Ratio",
             y="Edge Rank Threshold",
             title=title) +
        theme_ridges() +
        theme(axis.title.x=element_text(hjust=0.5),
              axis.title.y=element_text(hjust=0.25),
              plot.title=element_text(hjust=0.5))
}

run_int_ext_shuffle <- function(edgelist, geneset, rank_thresholds,
                                ntrials=1000, niter=10, epoch=250,
                                do_parallel=F) {

    ldply(rank_thresholds, function(rank_threshold) {

        cat("Running trials at rank threshold", rank_threshold, "\n")

        g <- edgelist %>% filter(Rank <= rank_threshold) %>%
            select(from=Gene.x, to=Gene.y) %>%
            graph_from_data_frame(directed=F,
                                  vertices=data_frame(Gene = unique(geneset$Gene)))

        graph_shuffle(g, geneset, "ratio", ntrials, niter,
                      epoch=epoch, do_parallel=do_parallel) %>%
            mutate(Rank = rank_threshold)
        #     group_by(Rank, Geneset) %>%
        #     summarise(IntExtShuffle = list(IntExtShuffle))

    }, .id=NA)
}


calc_int_ext_dens_ratio <-
    function(g, geneset_df, metric="enrichment",
             total_internal=NULL, total_external=NULL) {

        all_v <- names(V(g))


        # geneset_pairs <- geneset_df %>%
        #     inner_join(., ., by="Geneset") %>%
        #     filter(Gene.x < Gene.y)

        edgelist <-
            as_edgelist(g) %>%
            as.data.frame() %>%
            set_colnames(c("Gene.x", "Gene.y")) %>%
            mutate(Gene.x = as.character(Gene.x),
                   Gene.y = as.character(Gene.y)) %>%
            bind_rows(.,
                      .[,2:1] %>%
                          set_colnames(c("Gene.x", "Gene.y"))) %>%
            full_join(geneset_df, by=c("Gene.x" = "Gene"))

        int_ext <- edgelist %>%
            group_by(Geneset) %>%
            summarise(internal_vertices=intersect(Gene.x, all_v) %>% length,
                      external_vertices=length(all_v)-internal_vertices,
                      internal_true=sum(Gene.y %in% Gene.x, na.rm=T)/2,
                      external_true=sum(!Gene.y %in% Gene.x, na.rm=T)) %>%
            mutate(internal_total = choose(internal_vertices, 2),
                   external_total = internal_vertices * external_vertices,
                   internal_false = internal_total - internal_true,
                   external_false = external_total - external_true)



        if (metric == "ratio") {

            int_ext_metric <-
                ifelse(int_ext$internal_total == 0 |
                           int_ext$external_true == 0,
                       0,
                       (int_ext$internal_true / int_ext$internal_total) *
                           (int_ext$external_total / int_ext$external_true))

            # int_density <- int_ext$internal_true/int_ext$internal_total
            # ext_density <- int_ext$external_true/int_ext$external_total
            # int_ext_metric <- ifelse(ext_density==0, 0,
            #                          int_density/ext_density)

        } else if (metric == "enrichment") {

            int_ext_metric <- (int_ext$internal_true/int_ext$internal_false) *
                (int_ext$external_false/int_ext$external_true)
        } else {
            stop("Unknown metric")
        }
        return(int_ext_metric %>% set_names(int_ext$Geneset))

    }




graph_shuffle <- function(g, geneset_df, metric, ntrials, niter,
                          do_parallel=F, epoch=100) {


    # trials <- 1:ntrials
    trials_list <- split(1:ntrials, (0:(ntrials-1) %/% epoch))

    result_matrix <-
        matrix(0, nrow=ntrials, ncol=length(unique(geneset_df$Geneset))) %>%
        set_colnames(sort(unique(geneset_df$Geneset)))


    for (trials in trials_list) {
        if (trials[1] == 1) {
            cat("Beginning shuffle trials...\n")
        }

        t0 <- proc.time()

        result_matrix[trials,] <- laply(trials, function(i) {

            # return(data_frame(Trial=i, Test=rnorm(10)))
            # shuffle the graph
            g_shuffle <- rewire(g, keeping_degseq(loops=F, niter=niter*length(E(g))))


            # calc int/ext density ratio
            shuffled_int_ext <-
                calc_int_ext_dens_ratio(g_shuffle,
                                        geneset_df, metric=metric)

            # data_frame(Trial = i,
            #            Geneset=names(shuffled_int_ext),
            #            IntExtShuffle=shuffled_int_ext)

            return(shuffled_int_ext)

        }, .parallel=do_parallel)

        t1 <- proc.time()

        cat("Trial", tail(trials, 1), "in", t1[3]-t0[3], "seconds\n")
    }

    data_frame(Geneset = colnames(result_matrix),
               IntExtShuffle = alply(result_matrix, 2))
}
