
source("libraries.R")
source("upload_module.R")
source("volcano_module.R")
source("pca_module.R")   


ui <- navbarPage(
    
    title = "Seq2Viz",
    theme = bs_theme(bootswatch = "flatly"),
    
    # Welcome tab
    tabPanel(
        "Home",
        fluidPage(
            titlePanel("Welcome"),
            fluidRow(
                column(
                    width = 8,
                    p("This app allows you to explore your data through interactive plots."),
                    p("You can upload your files on the first page, then navigate to:"),
                    tags$ul(
                        tags$li("Volcano plot tab to visualize differential expression."),
                        tags$li("PCA plot tab to explore sample clustering.")
                    ),
                    p("Start by uploading your data on the Upload tab!")
                )
            )
        )
    ),
    
    tabPanel("Upload files",
             uploadUI("upload1")
    ),
    
    tabPanel("Volcano Plot",
             volcanoUI("volcano1")
    ),

    tabPanel("PCA Plot",
             pcaUI("pca1")
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
