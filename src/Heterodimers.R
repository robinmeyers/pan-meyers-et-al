library(ProjectTemplate)
load.project(override.config = list(cache_loading=F, munging=F))

#Heterodimers
load("./cache/rnai_dep_rank.RData")
load("./cache/avana_2017_dep_rank.RData")
load("./cache/coxpres_db.RData")
load("./cache/pdb_surfaces.RData")

out_dir <- file.path("./output/heterodimers", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)

rnai_ranked_df <-
    rnai_dep_rank %>%
    mat.to.df("Gene.x", "Gene.y", "Rank") %>%
    filter(Gene.x < Gene.y)

crispr_ranked_df <-
    avana_2017_dep_rank %>%
    mat.to.df("Gene.x", "Gene.y", "Rank") %>%
    filter(Gene.x < Gene.y)

coxpres_ranked_df <-
    coxpres_db %>%
    edgeweight_symmetric_rank %>%
    mat.to.df("Gene.x", "Gene.y", "Rank") %>%
    filter(Gene.x < Gene.y)


#PDB
surf_func <- pdb_surfaces %>%
    filter(Interaction > 500) %>%
    mutate(Max_Ratio = pmax(Interaction/Surface1, Interaction/Surface2)) %>%
    filter(Max_Ratio > 0.4) %>%
    inner_join(rnai_ranked_df, by = c("Gene.x", "Gene.y")) %>%
    inner_join(crispr_ranked_df, by = c("Gene.x", "Gene.y")) %>%
    inner_join(coxpres_ranked_df, by = c("Gene.x", "Gene.y")) %>%
    rename(RNAi = Rank.x,
           CRISPR = Rank.y,
           GE = Rank)

#scale by size
surf_func <- surf_func %>%
    mutate(RNAi = RNAi/nrow(rnai_dep_rank),
           CRISPR = CRISPR/nrow(avana_2017_dep_rank),
           GE = GE/nrow(coxpres_db))

#k-s test
ks.test(surf_func$RNAi, surf_func$CRISPR)
ks.test(surf_func$RNAi, surf_func$GE)
ks.test(surf_func$CRISPR, surf_func$GE)

data.frame(rank = seq(0, 1, 0.001)) %>%
    mutate(RNAi = ecdf(surf_func$RNAi)(rank),
           CRISPR = ecdf(surf_func$CRISPR)(rank),
           COXPRESdb = ecdf(surf_func$GE)(rank)) %>%
    gather(Dataset, ECDF, 2:4) %>%
    ggplot(aes(rank, ECDF, color = Dataset)) +
    geom_line() +
    labs(x = "Rank Correlation",
         y = "Frac. of Structural Heterodimer Recall") +
    scale_color_manual(values = c("#3B4CA0", "#F8766D", "#00BFC4", "#C77CFF")) +
    theme(legend.position = c(0,1),
          legend.justification = c(0,1))
    ggsave(file.path(out_dir, "heterodimers.pdf"),
           width = 4.5, height = 3)
