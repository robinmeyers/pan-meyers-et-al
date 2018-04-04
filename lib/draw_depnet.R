draw_depnet <- function(dep_corr, geneset_df, n_edges=5, seed=1234,
                        rank_method="maximal",
                        min_correlation=0,
                        min_connected_components=1,
                        allow_external_nodes=T, external_min_degree=2,
                        grow_smoothly=T,
                        add_labels=F,
                        internal_node_size=3, external_node_size=2,
                        plot_interactive=F,
                        plot_algorithm = "fr",
                        min_edge_width=.1,
                        ...) {

    set.seed(seed)

    geneset_df <- geneset_df %>%
        filter(Gene %in% rownames(dep_corr)) %>%
        mutate(Geneset = str_sub(Geneset, 1,20))

    if (length(unique(geneset_df$Geneset)) > 40) stop("Cannot handle more than 40 unique genesets")

    if ("Color" %in% colnames(geneset_df)) {
        node_colors <- geneset_df %>%
            distinct(Geneset, Color)  %$%
            {Color %>% set_names(Geneset)}
    } else {
        node_colors <- brewer_pal(palette="Paired")(10)[c(2,4,6,8,10,1,3,5,7,9)] %>%
            rep(4) %>% magrittr::extract(1:length(unique(geneset_df$Geneset)))
    }

    if ("Shape" %in% colnames(geneset_df)) {
        node_shapes <- geneset_df %>%
            distinct(Geneset, Shape)  %$%
            {Shape %>% set_names(Geneset)}
    } else {
        node_shapes <- rep(21:24, each=10) %>%
            magrittr::extract(1:length(unique(geneset_df$Geneset)))
    }

    diag(dep_corr) <- NA

    dep_corr[dep_corr <= min_correlation] <- NA


    if (allow_external_nodes) {

        edge_rank <- dep_corr
        edge_rank[!rownames(edge_rank) %in% geneset_df$Gene,] <- NA
        edge_rank[rownames(edge_rank) %in% geneset_df$Gene,] <-
            matrixStats::rowRanks(-edge_rank[rownames(edge_rank) %in% geneset_df$Gene,])



    } else {

        edge_rank <- dep_corr[rownames(dep_corr) %in% geneset_df$Gene,
                              rownames(dep_corr) %in% geneset_df$Gene]
        edge_rank <- matrixStats::rowRanks(-edge_rank) %>%
            set_rownames(rownames(edge_rank)) %>%
            set_colnames(colnames(edge_rank))


    }

    if (rank_method=="maximal") {
        edge_rank <- pmin(edge_rank, t(edge_rank), na.rm=T)
    }
    if (rank_method=="minimal") {
        edge_rank <- pmax(edge_rank, t(edge_rank), na.rm=T)
    }
    if (rank_method=="average") {
        edge_rank <- (edge_rank + t(edge_rank)) / 2
    }

    if (allow_external_nodes) {
        edge_rank <- edge_rank[rownames(edge_rank) %in% geneset_df$Gene,]
    }


    growing_graph <-
        make_empty_graph(directed=F) %>%
        add_vertices(length(rownames(edge_rank)), name=rownames(edge_rank))

    graph_coords <- matrix(rnorm(ncol(edge_rank)*2), nrow=ncol(edge_rank), ncol=2) %>%
        set_colnames(c("x", "y")) %>%
        set_rownames(colnames(edge_rank))

    for (i in 1:n_edges) {

        # add i-th best edge to each node

        edge_rank <-
            edge_rank[aaply(edge_rank, 1, function(x) !all(is.na(x))),]

        # edge_weights <- aaply(edge_matrix, 1, max, na.rm=T)

        next_closest_nodes <- which(edge_rank <= i, arr.ind=T)
        new_vertices <- setdiff(c(rownames(edge_rank)[next_closest_nodes[,1]],
                                  colnames(edge_rank)[next_closest_nodes[,2]]),
                                V(growing_graph)$name)
        new_edges_nonunique <- cbind(rownames(edge_rank)[next_closest_nodes[,1]],
                                     colnames(edge_rank)[next_closest_nodes[,2]])
        new_edges <- new_edges_nonunique[!duplicated(aaply(new_edges_nonunique, 1, sort)),]
        new_edge_weights <- 1/edge_rank[new_edges]

        # next_closest_nodes <- aaply(edge_rank, 1,
        #                            function(x) names(x)[which(x<=n_edges)[1]])
        # new_vertices <- setdiff(next_closest_node, V(growing_graph)$name)
        # new_edges <- cbind(names(next_closest_node), next_closest_node)
        # new_edge_weights <- 1/i


        edge_rank[new_edges_nonunique] <- NA

        growing_graph <- growing_graph %>%
            add_vertices(length(new_vertices), name=new_vertices) %>%
            add_edges(as.vector(t(new_edges)), weight=new_edge_weights)

        if (!grow_smoothly & i != n_edges) next


        growing_graph_tmp <- growing_graph %>%
            set_vertex_attr(name="Geneset",
                            value=data_frame(Gene=V(.)$name) %>%
                                left_join(geneset_df, by="Gene") %>%
                                distinct(Gene, .keep_all=T) %$%
                                as.character(Geneset)) %>%
            delete_vertices(which(is.na(V(.)$Geneset) & degree(.) < external_min_degree))

        connected_components <- components(growing_graph_tmp) %>%
        {.$csize[.$membership]}


        if (!all(connected_components < min_connected_components)) {
            growing_graph_tmp <- growing_graph_tmp %>%
                delete_vertices(which(connected_components < min_connected_components))
        }




        # delete_vertices(which(is.na(V(.)$Geneset) &
        #                           (adjacent_vertices(., V(.)) %>%
        #                                laply(length) %>%
        #                                is_less_than(external_min_degree))))




        if (plot_algorithm == "fr") {
            plotting_graph <- growing_graph_tmp %>%
                create_layout("igraph", algorithm='fr', coords=graph_coords[V(.)$name,],
                              ...)
        } else if (plot_algorithm == "kk") {
            plotting_graph <- growing_graph_tmp %>%
                create_layout("igraph", algorithm='kk', coords=graph_coords[V(.)$name,],
                              weights=1/V(.)$weight, ...)
        } else if (plot_algorithm == "graphopt") {
            plotting_graph <- growing_graph_tmp %>%
                create_layout("igraph", algorithm='graphopt', start=graph_coords[V(.)$name,],
                              ...)
        } else {
            stop("Error: don't recognize layout alogrithm")
        }


        graph_coords[as.character(plotting_graph$name),] <-
            as.matrix(plotting_graph[,c("x","y")])

    }

    growing_graph <- growing_graph %>%
        set_vertex_attr(name="Geneset",
                        value=data_frame(Gene=V(.)$name) %>%
                            left_join(geneset_df, by="Gene") %>%
                            distinct(Gene, .keep_all=T) %$%
                            as.character(Geneset) %>% str_sub(1,25))



    # growing_graph <- growing_graph %>%
    #     delete_vertices(which(is.na(V(.)$Geneset) &
    #                               (adjacent_vertices(., V(.)) %>%
    #                                    laply(length) %>%
    #                                    is_less_than(external_min_degree))))


    # growing_graph <- growing_graph %>%
    #     delete_vertices(which(adjacent_vertices(., V(.)) %>%
    #                               laply(length) %>%
    #                               is_less_than(internal_min_degree)))


    growing_graph <- growing_graph %>%
        delete_vertices(which(is.na(V(.)$Geneset) & degree(.) < external_min_degree))

    connected_components <- components(growing_graph) %>%
    {.$csize[.$membership]}

    if (!all(connected_components < min_connected_components)) {
        growing_graph <- growing_graph %>%
            delete_vertices(which(connected_components < min_connected_components))
    } else {
        return(list(graph=growing_graph, layout=plotting_graph, gg=NULL))
    }

    plotting_graph <- growing_graph %>%
        create_layout("manual",
                      node.positions=as.data.frame(graph_coords[V(.)$name,]))

    # plotting_graph_last <- plotting_graph

    # graph_coords[as.character(plotting_graph$name),] <-
    #     as.matrix(plotting_graph[,c("x","y")])

    max_width <- 1
    min_width <- max_width*max(min_edge_width, min(E(attr(plotting_graph, "graph"))$weight))

    gg <-
        # plotting_graph %>%
        ggraph(growing_graph, "manual",
               node.positions=as.data.frame(graph_coords[V(growing_graph)$name,])) +
        geom_edge_link0(aes(edge_width=ifelse(weight < min_edge_width, min_edge_width, weight)),
                        edge_colour="grey70") +
                        {if (plot_interactive) {
                            geom_point_interactive(aes(x=x, y=y, fill=Geneset, shape=Geneset,
                                                       size=as.character(!is.na(Geneset)),
                                                       tooltip=name, data_id=name),
                                                   alpha=1, stroke=0,
                                                   data=plotting_graph %>% arrange(!is.na(Geneset)))
                        } else {
                            geom_node_point(aes(fill=Geneset, shape=Geneset,
                                                size=as.character(!is.na(Geneset))),
                                            alpha=1, stroke=0,
                                            data=plotting_graph %>% arrange(!is.na(Geneset)))
                        }}

    if (add_labels) {
        label_df <- plotting_graph %>% filter(!is.na(Geneset))
        gg <- gg +
            geom_text_repel(aes(label=name, x=x, y=y),
                            data=label_df,
                            force=1,
                            size=2, segment.alpha=0,
                            box.padding=unit(0.1, "lines"),
                            # point.padding=unit(0.05, "lines"),
                            nudge_x=0.01*(label_df$x - mean(plotting_graph$x)) /
                                sqrt((label_df$x-mean(plotting_graph$x))^2 +
                                         (label_df$y-mean(plotting_graph$y))^2),
                            nudge_y=0.01*(label_df$y - mean(plotting_graph$y)) /
                                sqrt((label_df$x-mean(plotting_graph$x))^2 +
                                         (label_df$y-mean(plotting_graph$y))^2))
    }

    gg <- gg +
        scale_fill_manual(values=node_colors,
                          na.value="grey50") +
        # scale_color_brewer(type="qual", palette="Paired", na.value="grey50") +
        scale_shape_manual(values=node_shapes,
                           na.value=16) +
        # scale_edge_alpha_continuous(range=c(.05, 0.2), guide=F) +
        scale_edge_width_continuous(range=c(min_width, max_width), guide=F) +
        scale_size_manual(values=c("FALSE"=external_node_size,
                                   "TRUE"=internal_node_size), guide=F) +
        # labs(title=str_c(internal_node_size, "/", external_node_size)) +
        # labs(title=str_c(i,"-th best edge")) +
        # scale_colour_discrete() +
        # geom_node_text(aes(label = name, filter = degree > 150), color = 'white', size = 3) +
        ggforce::theme_no_axes() +
        guides(fill=guide_legend(ncol=1, override.aes=list(size=3), keyheight=unit(0.75, "lines"))) +
        theme(legend.title=element_blank(),
              legend.text=element_text(size=10),
              panel.border=element_blank(),
              legend.background=element_rect(color="grey70", size=0.5))

    # return(gg)

    return(list(graph=growing_graph, layout=plotting_graph, gg=gg))
}
