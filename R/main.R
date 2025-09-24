
source("libraries.R")
source("upload.R")
source("volcano.R")
source("pca.R")   
source("home.R")


ui <- navbarPage(
    
    title = "Seq2Viz",
    
    theme = bs_theme(
        preset = "cosmo",
        bg = "#F3F6FA",
        fg = "#2b5769",
        base_font = font_google("Jost"),
        # "navbar-bg" = "#2b5769",   
        "navbar-fg" = "#F3F6FA"
    ),



    tabPanel(
        "Home",
        homeUI("home1")
    ),
    
    
    tabPanel("Upload files",
             uploadUI("upload1")
    ),
    
    tabPanel("PCA Plot",
             pcaUI("pca")
    ),
    
    tabPanel("Volcano Plot",
             volcanoUI("volcano")
    ),

    tabPanel("Heatmap"
             # heatmapUI("heatmap")
    )
)

server <- function(input, output, session) {
    
    # returns a list of reactive datasets
    data_list <- callModule(uploadServer, "upload1")
    
    
    callModule(volcanoServer, "volcano", data_list$deseq)
    
    
    callModule(pcaServer, "pca", data_list$metadata, data_list$counts, data_list$deseq)
    
}


shinyApp(ui, server)
