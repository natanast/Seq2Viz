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
    ),
    
    tabPanel("Volcano Plot",
             sidebarLayout(
                 sidebarPanel(
                     numericInput("logfc_cutoff", "Log2FC cutoff", 1),
                     numericInput("padj_cutoff", "Adjusted p-value cutoff", 0.05)
                 ),
                 mainPanel(
                     plotOutput("volcano_plot", height = "600px")
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
    
    output$volcano_plot <- renderPlot({
        req(deseq_data())
        df <- deseq_data()
        
        req(all(c("log2FoldChange", "padj") %in% colnames(df)))
        
        df$Significant <- with(df,
                               ifelse(padj < input$padj_cutoff & abs(log2FoldChange) > input$logfc_cutoff,
                                      "Significant", "Not Significant"))
        
        ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
            geom_point(alpha = 0.6) +
            scale_color_manual(values = c("gray", "red")) +
            theme_minimal(base_size = 14) +
            labs(title = "Volcano Plot", x = "log2(Fold Change)", y = "-log10(Adjusted p-value)") +
            theme(legend.position = "top")
    })
}

shinyApp(ui, server)
