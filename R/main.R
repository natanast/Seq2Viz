

library(shiny)
library(data.table)


options(shiny.maxRequestSize = 50*1024^2)  # 50 MB limit


ui <- fluidPage(
    titlePanel("Seq2Viz: RNA-seq Analysis and Visualization"),
    
    sidebarLayout(
        sidebarPanel(
            h4("Load data"),
            fileInput("counts_file", "Upload counts.txt", 
                      accept = c(".txt", ".tsv")),
            helpText("Upload a counts file with gene names in the first column and sample names as column names.")
        ),
        
        mainPanel(
            h4("Preview of Uploaded Count Matrix"),
            tableOutput("counts_preview")
        )
    )
)

server <- function(input, output, session) {
    
    counts_data <- reactive({
        req(input$counts_file)
        fread(input$counts_file$datapath)
    })
    
    output$counts_preview <- renderTable({
        head(counts_data(), 10)
    })
    
}

shinyApp(ui, server)
