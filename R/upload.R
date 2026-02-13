
library(shiny)
library(data.table)
library(readxl)
library(tools)

uploadUI <- function(id) {
    ns <- NS(id)
    sidebarLayout(
        sidebarPanel(
            h4("1. Required for Analysis (DE Tab)"),
            fileInput(ns("metadata_file"), "Upload Sample Metadata (.xlsx, .csv)"),
            fileInput(ns("counts_file"), "Upload Raw Counts (.txt, .csv)"),
            
            hr(),
            
            h4("2. Optional (For Visualizations Only)"),
            p("Upload these if you already have results and want to skip the analysis step."),
            fileInput(ns("deseq_file"), "Upload DESeq2 Results (.xlsx, .csv)"),
            # --- ADDED THIS INPUT ---
            fileInput(ns("norm_counts_file"), "Upload Normalized Counts (.txt, .csv)"),
            
            style = "border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem; background-color: #fff;"
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("Metadata", tableOutput(ns("metadata_preview"))),
                tabPanel("Raw Counts", tableOutput(ns("counts_preview"))),
                tabPanel("DESeq Results", tableOutput(ns("deseq_preview"))),
                tabPanel("Norm Counts", tableOutput(ns("norm_preview")))
            ),
            style = "border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem; background-color: #fff;"
        )
    )
}

uploadServer <- function(id) {
    moduleServer(id, function(input, output, session) {
        
        # Helper to read files
        read_file <- function(file_input) {
            if (is.null(file_input)) return(NULL)
            ext <- tools::file_ext(file_input$name)
            tryCatch({
                if (ext %in% c("xlsx", "xls")) {
                    dt <- as.data.table(readxl::read_xlsx(file_input$datapath))
                } else {
                    # Auto-detect separator
                    dt <- fread(file_input$datapath)
                }
                return(dt)
            }, error = function(e) { return(NULL) })
        }
        
        # Reactive Data Loaders
        metadata_data <- reactive({ read_file(input$metadata_file) })
        counts_data   <- reactive({ read_file(input$counts_file) })
        deseq_data    <- reactive({ read_file(input$deseq_file) })
        # --- NEW READER ---
        norm_data     <- reactive({ read_file(input$norm_counts_file) })
        
        # Previews
        output$metadata_preview <- renderTable({ head(metadata_data(), 5) })
        output$counts_preview   <- renderTable({ head(counts_data(), 5) })
        output$deseq_preview    <- renderTable({ head(deseq_data(), 5) })
        output$norm_preview     <- renderTable({ head(norm_data(), 5) })
        
        # Return ALL data
        list(
            metadata = metadata_data,
            counts = counts_data,
            deseq = deseq_data,
            norm_counts = norm_data # <--- Return this new file
        )
    })
}

# library(shiny)
# library(data.table)
# library(readxl)
# library(tools)
# 
# uploadUI <- function(id) {
#     ns <- NS(id)
#     sidebarLayout(
#         sidebarPanel(
#             h4("1. Required for Analysis"),
#             fileInput(ns("metadata_file"), "Upload Sample Metadata (.xlsx, .csv)", 
#                       accept = c(".xlsx", ".csv", ".txt")),
#             
#             fileInput(ns("counts_file"), "Upload Raw Counts (.txt, .csv)", 
#                       accept = c(".csv", ".txt", ".tsv")),
#             
#             hr(),
#             
#             h4("2. Optional (If skipping DE)"),
#             fileInput(ns("deseq_file"), "Upload DESeq2 Results (.xlsx, .csv)", 
#                       accept = c(".xlsx", ".csv")),
#             
#             style = "border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem; background-color: #fff;"
#         ),
#         mainPanel(
#             tabsetPanel(
#                 tabPanel("Metadata Preview", 
#                          tableOutput(ns("metadata_preview"))
#                 ),
#                 tabPanel("Counts Preview", 
#                          tableOutput(ns("counts_preview"))
#                 ),
#                 tabPanel("DESeq2 Results Preview", 
#                          tableOutput(ns("deseq_preview"))
#                 )
#             ),
#             style = "border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem; background-color: #fff;"
#         )
#     )
# }
# 
# uploadServer <- function(id) {
#     moduleServer(id, function(input, output, session) {
#         
#         # --- Helper Function: Read Files Safely ---
#         read_file <- function(file_input) {
#             # If no file is uploaded, return NULL immediately
#             if (is.null(file_input)) return(NULL)
#             
#             ext <- tools::file_ext(file_input$name)
#             
#             tryCatch({
#                 if (ext %in% c("xlsx", "xls")) {
#                     dt <- as.data.table(readxl::read_xlsx(file_input$datapath))
#                 } else if (ext == "csv") {
#                     dt <- fread(file_input$datapath)
#                 } else if (ext %in% c("tsv", "txt")) {
#                     dt <- fread(file_input$datapath, sep = "\t") # Auto-detect separator usually works, but specific is safer
#                 } else {
#                     validate(paste("Unsupported file format:", file_input$name))
#                 }
#                 return(dt)
#             }, error = function(e) {
#                 showNotification(paste("Error reading file:", e$message), type = "error")
#                 return(NULL)
#             })
#         }
#         
#         # --- Reactive Data Loaders ---
#         
#         # 1. Metadata (Required usually)
#         metadata_data <- reactive({
#             read_file(input$metadata_file)
#         })
#         
#         # 2. Counts (Required for DE)
#         counts_data <- reactive({
#             read_file(input$counts_file)
#         })
#         
#         # 3. DESeq Results (Optional)
#         deseq_data <- reactive({
#             read_file(input$deseq_file)
#         })
#         
#         # --- Previews ---
#         output$metadata_preview <- renderTable({
#             req(metadata_data())
#             head(metadata_data(), 10)
#         })
#         
#         output$counts_preview <- renderTable({
#             req(counts_data())
#             # Show first 5 cols to keep it readable
#             df <- counts_data()
#             if(ncol(df) > 5) df[, 1:5, with=FALSE] else df
#             head(df, 10)
#         })
#         
#         output$deseq_preview <- renderTable({
#             req(deseq_data())
#             head(deseq_data(), 10)
#         })
#         
#         # --- Return Data to Main App ---
#         # We return the REACTIVES themselves, not the values
#         list(
#             metadata = metadata_data,
#             counts = counts_data,
#             deseq = deseq_data
#         )
#     })
# }


 
# uploadUI <- function(id) {
#     ns <- NS(id)
#     sidebarLayout(
#         sidebarPanel(
#             fileInput(ns("deseq_file"), "Upload DESeq2 results (.xlsx, .csv, .tsv)"),
#             fileInput(ns("metadata_file"), "Upload Sample Metadata (.xlsx, .csv, .tsv)"),
#             fileInput(ns("counts_file"), "Upload Counts (.txt, .csv, .tsv)"),
#             style = "border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem; background-color: #fff;"
#         ),
#         mainPanel(
#             h4("DESeq2 Results Preview"),
#             tableOutput(ns("deseq_preview")),
#             hr(),
#             h4("Sample Metadata Preview"),
#             tableOutput(ns("metadata_preview")),
#             hr(),
#             h4("Counts Preview"),
#             tableOutput(ns("counts_preview")),
#             style = "border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem; background-color: #fff;"
#         )
#     )
# }
# 
# 
# 
# uploadServer <- function(input, output, session) {
#     
#     # Helper function to read files
#     read_file <- function(file) {
#         req(file)
#         ext <- tools::file_ext(file$name)
#         if (ext %in% c("xlsx", "xls")) {
#             dt <- as.data.table(readxl::read_xlsx(file$datapath))
#         } else if (ext == "csv") {
#             dt <- fread(file$datapath)
#         } else if (ext %in% c("tsv", "txt")) {
#             dt <- fread(file$datapath, sep = "\t")
#         } else {
#             validate(paste("Unsupported file format:", file$name))
#         }
#         dt
#     }
#     
#     deseq_data <- reactive({ read_file(input$deseq_file) })
#     metadata_data <- reactive({ read_file(input$metadata_file) })
#     counts_data <- reactive({ read_file(input$counts_file) })
#     
#     output$deseq_preview <- renderTable({ head(deseq_data(), 10) })
#     output$metadata_preview <- renderTable({ head(metadata_data(), 10) })
#     output$counts_preview <- renderTable({ head(counts_data(), 10) })
#     
#     # Return a list of the three reactive
#     list(
#         deseq = deseq_data,
#         metadata = metadata_data,
#         counts = counts_data
#     )
# }
