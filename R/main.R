
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
    
    tabPanel("Volcano Plot",
             volcanoUI("volcano1")
    ),

    tabPanel("PCA Plot",
             pcaUI("pca1")
    ),
    
    tabPanel("Heatmap"
             # pcaUI("pca1")
    )
)

server <- function(input, output, session) {
    
    # Call Upload module: returns a list of reactive datasets
    data_list <- callModule(uploadServer, "upload1")
    
    # Pass DESeq2 data to volcano plot module
    callModule(volcanoServer, "volcano1", data_list$deseq)
    
    # Pass metadata and counts to PCA module
    callModule(pcaServer, "pca1", data_list$metadata, data_list$counts, data_list$deseq)
    
}


shinyApp(ui, server)
