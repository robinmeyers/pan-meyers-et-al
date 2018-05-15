
### Install missing packages

installed_pkgs <- installed.packages()

pkgs <-  c("ProjectTemplate",
		   "doMC",
		   "devtools",
		   "rfigshare",
		   "magrittr",
		   "readr",
		   "stringr",
		   "plyr",
		   "tibble",
		   "tidyr",
		   "dplyr",
		   "ggplot2",
		   "cowplot",
		   "gridExtra",
		   "scales",
		   "doMC",
		   "RColorBrewer",
		   "matrixStats",
		   "shiny",
		   "igraph",
		   "ggraph",
		   "ggiraph",
		   "ggrepel",
		   "ggridges",
		   "ggsignif",
		   "gplots",
		   "pheatmap"
		   )

if (length(setdiff(pkgs, installed_pkgs)) > 0) {
    install.packages(pkgs = setdiff(pkgs, installed_pkgs))
}

### Download data from figshare

library(rfigshare)

figshare_id <- 6005297

data_dir <- "./data/raw"

figshare_article <- fs_details(figshare_id)

for (figshare_file in figshare_article$files) {
    cat("downloading", figshare_file$name, "\n")
    file_path <- file.path(data_dir, figshare_file$name)
    download.file(figshare_file$download_url, file_path)
}


library(ProjectTemplate)
load.project()

