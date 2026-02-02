
source("libraries.R")
source("upload.R")
source("de.R")
source("volcano.R")
source("pca.R")   
source("home.R")
source("heatmap.R")


# Allow uploads up to 100 MB
options(shiny.maxRequestSize = 100*1024^2)


ui <- navbarPage(
    
    title = "Seq2Viz",
    
    theme = bs_theme(
        preset = "cosmo",
        bg = "#F3F6FA",
        fg = "#2b5769",
        base_font = font_google("Jost"),
        "navbar-fg" = "#F3F6FA"
    ),


    tabPanel(
        "Home",
        homeUI("home1")
    ),
    
    tabPanel(
        "Upload files",
         uploadUI("upload1")
    ),
    
    tabPanel(
        "Deferential exprassion",
        deUI("de")
    ),
    
    tabPanel(
        "PCA Plot",
         pcaUI("pca")
    ),
    
    tabPanel(
        "Volcano Plot",
         volcanoUI("volcano")
    ),

    tabPanel(
        "Heatmap",
         heatmapUI("heat1")
    )
)



server <- function(input, output, session) {
    
    data_list <- callModule(uploadServer, "upload1")
    
    # callModule(deserver, "de")
    
    callModule(volcanoServer, "volcano", data_list$deseq)
    
    callModule(pcaServer, "pca", data_list$metadata, data_list$counts, data_list$deseq)
    
    heatmapServer(
        "heat1",
        meta_data   = data_list$metadata,
        counts_data = data_list$counts,
        deseq_data  = data_list$deseq
    )
    
    
}


shinyApp(ui, server)
