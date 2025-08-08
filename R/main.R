
source("libraries.R")
source("upload_module.R")
source("volcano_module.R")
source("pca_module.R")   # <-- Add this line


ui <- navbarPage(
    title = "Seq2Viz",
    theme = bs_theme(bootswatch = "flatly"),
    
    tabPanel("Upload DESeq2 Results",
             uploadUI("upload1")
    ),
    
    tabPanel("Volcano Plot",
             volcanoUI("volcano1")
    )
    # 
    # tabPanel("PCA Plot",      
    #          pca_ui("pca1")
    # )
)

server <- function(input, output, session) {
    # Call Upload module: returns a list of reactive datasets
    data_list <- callModule(uploadServer, "upload1")
    
    # data_list$deseq is reactive DESeq2 results
    # data_list$metadata is reactive metadata
    # data_list$counts is reactive counts
    
    # Pass DESeq2 data to volcano plot module
    callModule(volcanoServer, "volcano1", data_list$deseq)
    
    # You can now add more modules and pass data_list$metadata and data_list$counts as needed
}


shinyApp(ui, server)
