
cache("corum", {
    read_tsv("./data/raw/corum_human_core_complexes.tsv",
             col_types = "icc") %>%
        distinct(Complex, Gene, .keep_all = T) %>%
        add_count(Complex)
})

cache("corum_list", {
    split(corum$Gene, corum$Complex)
}, depends="corum")

# cache("c2_kegg_list", {
#     read.gmt("./data/raw/c2.cp.kegg.v6.0.symbols.gmt", as.df=F)
# })
#
# cache("c2_reactome_list", {
#     read.gmt("./data/raw/c2.cp.reactome.v6.0.symbols.gmt", as.df=F)
# })

cache("humap_list", {
    scan("./data/raw/humap_clusters.txt", what="", sep="\n")  %>%
        str_split(pattern = "\t") %>%
        set_names(1:length(.))
})

cache("baf_genes", {
    c("SMARCA4", "SMARCB1", "SMARCE1", "ARID1A", "SMARCC1", "SMARCC2",
      "SMARCD1", "SMARCD2", "BRD9", "GLTSCR1",
      "ARID2", "PBRM1", "BRD7", "PHF10",
      "SMARCA2", "ARID1B")
})

# cache("hart_2015_essentials", {
#     read.csv("./data/raw/hart_2015_essentials.csv", stringsAsFactors=FALSE)
# })
#
#
# cache("hart_nonessentials", {
#     hart_ctrls <- read_tsv("./data/raw/hart2014_seed_core_essentials.txt")
#     ne_genes <- hart_ctrls$`Nonessential Genes (NE)` %>% str_subset("[A-Z0-9]+")
# })
#
# cache("hart_coreessentials", {
#     hart_ctrls <- read_tsv("./data/raw/hart2014_seed_core_essentials.txt")
#     cce_genes <- hart_ctrls$`ConstitutiveCoreEssentials(CCE)` %>% str_subset("[A-Z0-9]+")
# })
