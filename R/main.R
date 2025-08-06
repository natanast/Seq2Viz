

library(shiny)
library(data.table)
library(readxl)
library(stringr)
library(ggplot2)
library(ggrepel)
library(colorspace)
library(shadowtext)
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
                     plotOutput("volcano_plot", height = "700px")
                 )
                 
             )
             
    )
    
)

server <- function(input, output, session) {
    
    deseq_data <- reactive({
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
        head(deseq_data(), 10)
    })
    
    output$volcano_plot <- renderPlot({
        
        df <- deseq_data()
        
        req(all(c("log2FoldChange", "padj", "GeneID") %in% colnames(df)))
        
        logfc_cutoff <- input$logfc_cutoff
        padj_cutoff <- input$padj_cutoff
        
        df <- df[!is.na(padj)]
        df$y <- -log10(df$padj)
        
        # Annotation
        df[, ann := fifelse(
            padj > padj_cutoff, "Not significant",
            fifelse(log2FoldChange > 0, "Up regulated", "Down regulated")
        )]
        
        df$ann = fifelse(
            df$padj <= 0.05 & df$log2FoldChange > -1 & df$log2FoldChange < 1, 
            paste0(df$ann, " (low)"),
            df$ann
        )
        
        # Top genes for labeling
        df2 <- df[padj <= padj_cutoff & abs(log2FoldChange) >= logfc_cutoff]
        df2 <- df2[order(abs(log2FoldChange), decreasing = TRUE)]
        df2 <- df2[, by = ann, head(.SD, 10) ]
        
        
        # Plot
        ggplot(df, aes(x = log2FoldChange, y = y)) +
            
            geom_point(aes(fill = ann),
                       shape = 21, stroke = NA, size = 2, alpha = 0.5) +
            
            geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff),
                       linewidth = 0.3, linetype = "dashed") +
            
            geom_hline(yintercept = -log10(padj_cutoff),
                       linewidth = 0.3, linetype = "dashed") +
            
            geom_point(data = df2, aes(x = log2FoldChange, y = y, fill = ann),
                       shape = 21, stroke = 0.15, size = 2, color = "white") +
            
            geom_text_repel(
                data = df2,
                aes(label = GeneID),
                max.overlaps = Inf,
                fontface = "bold",
                size = 2.5,
                bg.color = "white",
                bg.r = 0.05
            ) +
            
            scale_fill_manual(
                values = c(
                    "Up regulated" = "#990000",
                    "Up regulated (low)" = lighten("#990000", .5),
                    "Down regulated" = "#004d99",
                    "Down regulated (low)" = lighten("#004d99", .5),
                    "Not significant" = "grey"
                ),
                breaks = c("Up regulated", "Not significant", "Down regulated"),
                guide = guide_legend(override.aes = list(size = 3, alpha = 1))
            ) +
            
            scale_x_continuous(
                breaks = c(-5, -2.5, -1, 0, 1, 2.5, 5),
                trans = scales::pseudo_log_trans()
            ) +
            
            scale_y_continuous(
                expand = c(0, 0),
                breaks = c(2, 5, 10, 20, 30, 40),
                trans = scales::pseudo_log_trans()
            ) +
            
            coord_cartesian(clip = "off") +
            
            theme_minimal() +
            
            theme(
                legend.title = element_blank(),
                legend.position = "bottom",
                axis.line = element_line(linewidth = .3, color = "black"),
                axis.ticks = element_line(linewidth = .3, color = "black"),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_line(linewidth = .3, linetype = "dashed", color = "grey85"),
                plot.margin = margin(20, 20, 20, 20)
            ) +
            
            labs(x = "log2(Fold Change)", y = "-log10(padj)")
    })
    
}

shinyApp(ui, server)
