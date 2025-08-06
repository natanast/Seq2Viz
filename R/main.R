library(shiny)
library(data.table)
library(readxl)      # <-- NEW: For Excel
library(ggplot2)
library(ggrepel)
library(bslib)

ui <- navbarPage(
    title = "Seq2Viz",
    theme = bs_theme(bootswatch = "flatly"),
    
    tabPanel("Upload DESeq2 Results",
             sidebarLayout(
                 sidebarPanel(
                     fileInput("deseq_file", "Upload DESeq2 results (.xlsx, .csv, .tsv)")
                 ),
                 mainPanel(
                     tableOutput("deseq_preview")
                 )
             )
    )

)

server <- function(input, output, session) {
    
    deseq_data <- reactive({
        req(input$deseq_file)
        ext <- tools::file_ext(input$deseq_file$name)
        
        if (ext %in% c("xlsx", "xls")) {
            read_excel(input$deseq_file$datapath)
        } else if (ext == "csv") {
            fread(input$deseq_file$datapath)
        } else if (ext %in% c("tsv", "txt")) {
            fread(input$deseq_file$datapath, sep = "\t")
        } else {
            validate("Unsupported file format. Please upload .xlsx, .csv, or .tsv")
        }
    })
    
    output$deseq_preview <- renderTable({
        head(deseq_data(), 10)
    })
    

}

shinyApp(ui, server)
