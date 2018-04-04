
library(shiny)
library(igraph)
library(ggraph)
library(ggiraph)
library(ggrepel)

library(ProjectTemplate)
load.project(override.config = list(cache_loading=F,
                                    munging=F))

load("./cache/avana_dep_corr.RData")
load("./cache/rnai_dep_corr.RData")
load("./cache/coxpres_db.RData")
load("./cache/corum_list.RData")


custom_genesets <- list(BAF = c("SMARCA4", "SMARCA2", "ARID1A", "SMARCB1", "SMARCE1",
                                "DPF2", "SS18", "SMARCD1", "BRD9",
                                "PBRM1", "ARID2", "BRD7"),
                        PRC2 = c("SUZ12", "EZH2", "EED", "RBBP7"),
                        `RNA Pol II` = c("POLR2A", "POLR2B", "POLR2C", "POLR2D",
                                         "POLR2E", "POLR2F", "POLR2G", "POLR2H",
                                         "POLR2I", "POLR2J", "POLR2K", "POLR2L"),
                        Mediator = c("MED6", "MED8", "MED11", "MED17", "MED18", "MED20", "MED22",
                                     "MED27", "MED28", "MED29", "MED30", "MED1", "MED4", "MED7", "MED9",
                                     "MED10", "MED19", "MED21", "MED31", "MED26", "MED14", "MED15",
                                     "MED16", "MED23", "MED24", "MED25", "MED12", "MED12L", "MED13",
                                     "MED13L", "CDK8", "CDK19", "CCNC"),
                        Cohesin = c("RAD21", "REC8", "SMC1A", "SMC1B", "SMC3", "STAG1", "STAG2",
                                    "STAG3"),
                        ANAPC = c("ANAPC1", "ANAPC2", "ANAPC4", "ANAPC5",
                                  "ANAPC7", "ANAPC10", "ANAPC11", "ANAPC13", "ANAPC15",
                                  "CDC16", "CDC23", "CDC26", "CDC27"),
                        Histones = c("HIST2H3D", "HIST1H2AK", "HIST1H2AD", "HIST1H2AC",
                                     "HIST1H2BK", "HIST1H2BD", "HIST2H2BA", "HIST1H2BC",
                                     "HIST2H4B", "HIST2H4A", "HIST1H4J", "HIST1H4K",
                                     "HIST1H1E", "HIST1H2AH", "HIST1H3H", "HIST1H4I",
                                     "HIST1H2AE", "HIST1H2BG", "HIST1H4D"))

complex_geneset <- c(custom_genesets, corum_list) %>%
    data_frame(Geneset = names(.), Gene = .) %>%
    unnest(Gene)


# Define UI for application that draws a histogram
ui <- fluidPage(

   # Application title
   titlePanel("Protein Complex Co-Dependency Network"),

   ggiraphOutput("complex_depnet"),

   wellPanel(fluidRow(
       column(3,
              radioButtons("dataset", "", choices=c("RNAi", "CRISPR", "COXPRESdb"), inline=T),
              selectizeInput("complexes",
                          "Complexes:",
                          choices=unique(complex_geneset$Geneset),
                          selected=names(custom_genesets),
                          multiple=T),
              fileInput("geneset_file", "Upload Custom Set:")),
       column(3,
           numericInput("n_edges", "Number of edges:", value = 5,
                        min = 1, max = 50, step=1),
           numericInput("min_correlation", "Minimum Correlation:",
                          value=0.2, min=0, max=0.6, step=0.01),
           numericInput("min_connected_components", "Min Connected:", value=1,
                        min=1, step=1),
           checkboxInput("external_nodes", "External Nodes?", value=T),
           numericInput("external_min_degree", "External Node Degree:", value=2,
                        min=1, max=20, step=1)),
       column(2,
              numericInput("niter", "Iterations:", value=100,
                           min=1, max=100000, step=100),
              numericInput("start.temp", "Start temp:", value=1,
                           min=0.01, max=10, step=0.1),
              numericInput("internal_node_size", "Int Node Size:", value=3,
                           min=0, max=5, step=1),
              numericInput("external_node_size", "Ext Node Size:", value=2,
                           min=0, max=5, step=1)),
       column(2,
              numericInput("seed", "Seed:", value=1234,
                           min=0, max=9999, step=1),
              checkboxInput("add_label", "Gene Labels?", value=F),
              numericInput("figwidth", "Width:", value=6,
                           min=1, max=100, step=1),
              numericInput("figheight", "Height:", value=4,
                           min=1, max=100, step=1),
              downloadButton("save_pdf", "Save PDF"),
              downloadButton("save_png", "Save PNG"))
   ))


)

# Define server logic required to draw a histogram
server <- function(input, output) {


    depnetPlot <- reactive({

        if (input$dataset == "RNAi") {
            dep_corr <- rnai_dep_corr
        }
        else if (input$dataset == "CRISPR") {
            dep_corr <- avana_dep_corr
        } else {
            dep_corr <- coxpres_db
        }


        if (is.null(input$geneset_file)) {
            geneset_df <- complex_geneset %>%
                filter(Geneset %in% input$complexes)
        } else {
            geneset_df <- read_tsv(input$geneset_file$datapath)
        }

        depnet <- draw_depnet(dep_corr, geneset_df,
                              n_edges=ifelse(is.na(input$n_edges), 0, input$n_edges),
                              min_correlation=input$min_correlation,
                              grow_smoothly=T,
                              min_connected_components=input$min_connected_components,
                              allow_external_nodes=input$external_nodes,
                              external_min_degree=input$external_min_degree,
                              plot_interactive=T,
                              add_label=input$add_label,
                              internal_node_size=input$internal_node_size,
                              external_node_size=input$external_node_size,
                              niter=input$niter, start.temp=input$start.temp,
                              seed=input$seed)
        depnet_gg <- depnet$gg
    })


   output$complex_depnet <- renderggiraph({

      ggiraph(print(depnetPlot()), hover_css="fill:yellow",
              width_svg=6, height_svg=4, selection_type = "none")

   })

   output$save_pdf <- downloadHandler(
       filename = function() {
           'depnet.pdf'
       },
       content = function(file) {
           ggsave(file, depnetPlot(), width=input$figwidth, height=input$figheight)
       },
       contentType="application/pdf"
   )
   output$save_png <- downloadHandler(
       filename = function() {
           'depnet.png'
       },
       content = function(file) {
           ggsave(file, depnetPlot(), width=input$figwidth, height=input$figheight)
       },
       contentType="image/png"
   )

}

# Run the application
app <- shinyApp(ui = ui, server = server)
