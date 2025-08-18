
source("libraries.R")
source("upload_module.R")
source("volcano_module.R")
source("pca_module.R")   


ui <- navbarPage(
    
    title = "Seq2Viz",
    theme = bs_theme(bootswatch = "flatly"),
    
    tabPanel("Home",
             uploadUI("upload1")
    ),
    
    tabPanel("Upload files",
             uploadUI("upload1")
    ),
    
    tabPanel("Volcano Plot",
             volcanoUI("volcano1")
    ),

    tabPanel("PCA Plot",
             pcaUi("pca1")
    )
)

server <- function(input, output, session) {
    
    # Call Upload module: returns a list of reactive datasets
    data_list <- callModule(uploadServer, "upload1")
    
    # Pass DESeq2 data to volcano plot module
    callModule(volcanoServer, "volcano1", data_list$deseq)
    
    # Pass metadata and counts to PCA module
    callModule(pcaServer, "pca1", data_list$metadata, data_list$counts)
    
}


shinyApp(ui, server)
