library(ProjectTemplate)
load.project(override.config = list(cache_loading=F, munging=F))

load("./cache/avana_2017_dep.RData")
load("./cache/rnai_dep.RData")

out_dir <- file.path("./output/cell_line_tissue", Sys.Date())
dir.create(out_dir, recursive=T, showWarnings=F)


lineage_barplot <- function(dataset) {
  tibble(Samples = colnames(dataset)) %>% 
    mutate(Cell_Line = word(Samples, start = 1, end = 1, sep = "_"),
           Tissue = word(Samples, start = 2, end = -1, sep = "_")) %>% 
    group_by(Tissue) %>% 
    dplyr::count() %>% 
    ggplot(aes(Tissue, n)) +
    geom_col() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Lineage", y = "Number of Lines")
}

lineage_barplot(avana_2017_dep) +
  ggsave(file.path(out_dir, "crispr_tissue.pdf"),
         width = 8, height = 7)

lineage_barplot(rnai_dep) +
  ggsave(file.path(out_dir, "rnai_tissue.pdf"),
         width = 8, height = 7)
