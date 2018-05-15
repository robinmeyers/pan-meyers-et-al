library(ProjectTemplate)
load.project(override.config = list(cache_loading=F, munging=F))

library(ggjoy)
library(matrixStats)
library(pheatmap)

out_dir <- file.path("./output/fitness_landscapes", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)

load("./cache/avana_2017_dep.RData")

generate_landscape <- function(geneset, show_genes = TRUE, scale_input = 1) {
    #Data loading
    dat <- avana_2017_dep[intersect(rownames(avana_2017_dep),geneset),]

    #Data munging
    #1. Switch essentiality to positive
    #2. Scaling between 1 and 0
    #3. Soft thresholding (take to third power)

    dat <- -dat
    dat <- dat-rowMins(dat, na.rm = TRUE)
    dat <- dat/rowMaxs(dat, na.rm = TRUE)

    dat <- dat ^ 3

    ## Row clustering
    hr <- hclust(as.dist(1-cor(t(dat), method="pearson", use = "pairwise")),
                 method="complete")

    ## Column clustering
    hc <- hclust(as.dist(1-cor(dat, method="pearson", use = "pairwise")), method="complete")

    #Plot
    joy <- dat[rev(hr$labels[hr$order]), hc$labels[hc$order]] %>%
        mat.to.df(row.name = "Gene", col.name = "Cell_Line", dat.name = "Dep") %>%
        mutate(Gene = factor(Gene, levels = hr$labels[hr$order]),
               Cell_Line = factor(Cell_Line, levels = hc$labels[hc$order])) %>%
        ggplot(aes(x = Cell_Line, y = Gene, height = Dep, group = Gene)) +
        geom_joy(stat = "identity", scale = scale_input, size = 0) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank()) +
        labs(x = "Cell lines")
    if(show_genes)
        return(joy)
    else
        return(joy + theme(axis.text.y = element_blank(),
                           axis.ticks.y = element_blank()))
}

#Figure 1
fig1_geneset <- c("ELP6","ELP3","ELP2","ELP5","ELP4","EIF3J","ARPC3","ACTR3","ARPC4","ARPC2","ACTR2","EIF3I","EIF3D","EIF3A","EIF3L","EIF3K","NDUFA7","SDHB","SDHA","SDHC","SDHD","NDUFA3","NDUFA13","NDUFA8","NDUFA10","NDUFA1","NDUFA2","NDUFA11","NDUFA6","NDUFA9")
generate_landscape(fig1_geneset)
ggsave(file.path(out_dir, "overview_landscape.pdf"),
                 width = 5, height = 7)

#Figure 3
generate_landscape(c("MED24", "MED23", "MED16", "MED9", "MED1", "MED8", "MED18", "MED20"))
ggsave(file.path(out_dir, "med_crispr_landscape.pdf"),
       width = 5, height = 4)

#Others
#mSWI/SNF
c("SMARCA4", "SMARCB1", "SMARCE1", "ARID1A", "SMARCC1",
  "SMARCD1", "BRD9", "GLTSCR1",
  "ARID2", "PBRM1", "BRD7", "DPF2") %>%
    generate_landscape()

#STAGA/ATAC
read_tsv("./data/manual/staga_atac.txt") %>%
    pull(Gene) %>%
    generate_landscape()

#Mito
read_tsv("./data/manual/mito_respiratory_chain.txt") %>%
    pull(Gene) %>%
    generate_landscape()
