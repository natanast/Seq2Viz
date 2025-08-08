
source("libraries.R")

uploadUI <- function(id) {
    ns <- NS(id)
    sidebarLayout(
        sidebarPanel(
            fileInput(ns("deseq_file"), "Upload DESeq2 results (.xlsx, .csv, .tsv)")
        ),
        mainPanel(
            tableOutput(ns("deseq_preview"))
        )
    )
}

uploadServer <- function(input, output, session) {
    data <- reactive({
        req(input$deseq_file)
        ext <- tools::file_ext(input$deseq_file$name)
        
        if (ext %in% c("xlsx", "xls")) {
            as.data.table(readxl::read_xlsx(input$deseq_file$datapath))
        } else if (ext == "csv") {
            fread(input$deseq_file$datapath)
        } else if (ext %in% c("tsv", "txt")) {
            fread(input$deseq_file$datapath, sep = "\t")
        } else {
            validate("Unsupported file format. Please upload .xlsx, .csv, or .tsv")
        }
    })
    
    output$deseq_preview <- renderTable({
        head(data(), 10)
    })
    
    return(data)
}
