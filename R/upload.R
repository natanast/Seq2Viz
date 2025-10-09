
uploadUI <- function(id) {
    ns <- NS(id)
    sidebarLayout(
        sidebarPanel(
            fileInput(ns("deseq_file"), "Upload DESeq2 results (.xlsx, .csv, .tsv)"),
            fileInput(ns("metadata_file"), "Upload Sample Metadata (.xlsx, .csv, .tsv)"),
            fileInput(ns("counts_file"), "Upload Counts (.txt, .csv, .tsv)"),
            style = "border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem; background-color: #fff;"
        ),
        mainPanel(
            h4("DESeq2 Results Preview"),
            tableOutput(ns("deseq_preview")),
            hr(),
            h4("Sample Metadata Preview"),
            tableOutput(ns("metadata_preview")),
            hr(),
            h4("Counts Preview"),
            tableOutput(ns("counts_preview")),
            style = "border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem; background-color: #fff;"
        )
    )
}



uploadServer <- function(input, output, session) {
    
    # Helper function to read files
    read_file <- function(file) {
        req(file)
        ext <- tools::file_ext(file$name)
        if (ext %in% c("xlsx", "xls")) {
            dt <- as.data.table(readxl::read_xlsx(file$datapath))
        } else if (ext == "csv") {
            dt <- fread(file$datapath)
        } else if (ext %in% c("tsv", "txt")) {
            dt <- fread(file$datapath, sep = "\t")
        } else {
            validate(paste("Unsupported file format:", file$name))
        }
        dt
    }
    
    deseq_data <- reactive({ read_file(input$deseq_file) })
    metadata_data <- reactive({ read_file(input$metadata_file) })
    counts_data <- reactive({ read_file(input$counts_file) })
    
    output$deseq_preview <- renderTable({ head(deseq_data(), 10) })
    output$metadata_preview <- renderTable({ head(metadata_data(), 10) })
    output$counts_preview <- renderTable({ head(counts_data(), 10) })
    
    # Return a list of the three reactive
    list(
        deseq = deseq_data,
        metadata = metadata_data,
        counts = counts_data
    )
}
