
source("libraries.R")
source("upload_module.R")
source("volcano_module.R")

ui <- navbarPage(
    title = "Seq2Viz",
    theme = bs_theme(bootswatch = "flatly"),
    
    tabPanel("Upload DESeq2 Results",
             uploadUI("upload1")
    ),
    
    tabPanel("Volcano Plot",
             volcanoUI("volcano1")
    )
)

server <- function(input, output, session) {
    # Call Upload module: returns reactive data.table with DESeq2 data
    deseq_data <- callModule(uploadServer, "upload1")
    
    # Pass reactive data to Volcano plot module
    callModule(volcanoServer, "volcano1", deseq_data)
}

shinyApp(ui, server)
