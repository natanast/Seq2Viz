
library(shiny)
library(ggplot2)
library(data.table)
library(ggrepel)
library(colorspace)
library(shinycssloaders)

volcanoUI <- function(id) {
    ns <- NS(id)
    sidebarLayout(
        sidebarPanel(
            h4("Volcano Settings"),
            numericInput(ns("logfc_cutoff"), "Log2FC cutoff", 1, min = 0, step = 0.5),
            selectInput(ns("p_metric"), "Significance metric:",
                        choices = c("Adjusted p-value (padj)" = "padj",
                                    "p-value" = "pvalue")),
            numericInput(ns("p_cutoff"), "Significance Cutoff", 0.05, min = 0, max = 1, step = 0.01),
            hr(),
            downloadButton(ns("download_plot"), "Download Plot")
        ),
        mainPanel(
            # Added error message output for debugging
            verbatimTextOutput(ns("error_msg")),
            withSpinner(plotOutput(ns("volcano_plot"), height = "700px"))
        )
    )
}

volcanoServer <- function(id, deseq_data) {
    moduleServer(id, function(input, output, session) {
        
        # 1. Prepare Data (The Smart Renamer)
        clean_data <- reactive({
            req(deseq_data())
            
            # Create a copy and ensure it is a data.table
            df <- copy(deseq_data())
            if (!is.data.table(df)) setDT(df)
            
            # --- FIX 1: Normalize ID Column to "Geneid" ---
            # Your plot expects "Geneid". We find whatever exists and rename it to "Geneid".
            if (!"Geneid" %in% colnames(df)) {
                candidates <- c("gene_name", "gene_id", "id", "ID", "Gene", "Symbol", "row", "Row.names")
                found <- intersect(colnames(df), candidates)
                if (length(found) > 0) {
                    setnames(df, found[1], "Geneid")
                } else {
                    # Fallback: Use rownames if no ID column found
                    df$Geneid <- rownames(df)
                }
            }
            
            # --- FIX 2: Normalize LogFC Column ---
            if (!"log2FoldChange" %in% colnames(df)) {
                if ("logFC" %in% colnames(df)) setnames(df, "logFC", "log2FoldChange")
            }
            
            # --- FIX 3: Check P-value Columns ---
            # If user wants padj but file has 'adj.P.Val' (common in other tools), fix it
            if (!"padj" %in% colnames(df) && "adj.P.Val" %in% colnames(df)) {
                setnames(df, "adj.P.Val", "padj")
            }
            if (!"pvalue" %in% colnames(df) && "P.Value" %in% colnames(df)) {
                setnames(df, "P.Value", "pvalue")
            }
            
            df
        })
        
        # 2. Plotting Logic
        volcano_plot <- reactive({
            df <- clean_data()
            
            # Debugging Checks
            req(df)
            if (!"log2FoldChange" %in% colnames(df)) return(NULL)
            
            p_metric <- input$p_metric
            if (!p_metric %in% colnames(df)) {
                # If selected metric is missing, try to fallback to the other one
                if (p_metric == "padj" && "pvalue" %in% colnames(df)) {
                    p_metric <- "pvalue"
                    showNotification("Column 'padj' not found. Using 'pvalue' instead.", type = "warning")
                } else {
                    return(NULL)
                }
            }
            
            logfc_cutoff <- input$logfc_cutoff
            p_cutoff <- input$p_cutoff
            
            # Filter out NAs
            df <- df[!is.na(get(p_metric)) & !is.na(log2FoldChange)]
            
            # Calculate Y axis
            df$y <- -log10(df[[p_metric]])
            
            # Annotate
            df[, ann := fifelse(
                get(p_metric) > p_cutoff, "Not significant",
                fifelse(log2FoldChange > 0, "Up regulated", "Down regulated")
            )]
            
            # Add "low" annotation logic
            df$ann <- fifelse(
                df[[p_metric]] <= p_cutoff & abs(df$log2FoldChange) < logfc_cutoff & df$log2FoldChange != 0, 
                paste0(df$ann, " (low)"),
                df$ann
            )
            
            # Handle the case where "low" creates categories that might not have colors defined
            # Simplified logic: If it passes P-value but fails LogFC, mark as Not Significant for cleaner plot?
            # Or keep your specific logic. Let's stick to yours but ensure colors exist.
            
            # Subset for top genes (labels)
            # We explicitly require Geneid here, which we ensured exists in Step 1
            df2 <- df[get(p_metric) <= p_cutoff & abs(log2FoldChange) >= logfc_cutoff]
            df2 <- df2[order(abs(log2FoldChange), decreasing = TRUE)]
            
            # Safe head() in case df2 is empty
            if(nrow(df2) > 0) {
                df2 <- df2[, head(.SD, 10), by = ann]
            }
            
            # Plot
            p <- ggplot(df, aes(x = log2FoldChange, y = y)) +
                geom_point(aes(fill = ann), shape = 21, stroke = NA, size = 2, alpha = 0.5) +
                geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linewidth = 0.3, linetype = "dashed") +
                geom_hline(yintercept = -log10(p_cutoff), linewidth = 0.3, linetype = "dashed") +
                
                scale_fill_manual(
                    values = c(
                        "Up regulated" = "#990000",
                        "Up regulated (low)" = colorspace::lighten("#990000", .5),
                        "Down regulated" = "#004d99",
                        "Down regulated (low)" = colorspace::lighten("#004d99", .5),
                        "Not significant" = "grey"
                    ),
                    # Ensure all possible keys are covered or let ggplot drop unused ones
                    na.value = "grey"
                ) +
                scale_x_continuous(breaks = c(-5, -2.5, -1, 0, 1, 2.5, 5),
                                   trans = scales::pseudo_log_trans()) +
                scale_y_continuous(expand = c(0, 0),
                                   breaks = c(2, 5, 10, 20, 30, 40),
                                   trans = scales::pseudo_log_trans()) +
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
                labs(x = "log2(Fold Change)", y = paste0("-log10(", p_metric, ")"))
            
            # Add labels if df2 has data
            if (nrow(df2) > 0) {
                p <- p + 
                    geom_point(data = df2, aes(x = log2FoldChange, y = y, fill = ann), 
                               shape = 21, stroke = 0.15, size = 2, color = "white") +
                    geom_text_repel(data = df2, aes(label = Geneid), max.overlaps = Inf,
                                    fontface = "bold", size = 2.5, bg.color = "white", bg.r = 0.05)
            }
            
            p
        })
        
        # 3. Error Handling Output
        output$error_msg <- renderText({
            df <- clean_data()
            if (is.null(df)) return("Waiting for data...")
            
            # Check for critical columns
            missing <- c()
            if (!"log2FoldChange" %in% colnames(df)) missing <- c(missing, "log2FoldChange")
            if (!input$p_metric %in% colnames(df)) missing <- c(missing, input$p_metric)
            
            if (length(missing) > 0) {
                return(paste("Error: The following columns are missing from your data:", paste(missing, collapse=", ")))
            }
            return("")
        })
        
        output$volcano_plot <- renderPlot({ volcano_plot() })
        
        output$download_plot <- downloadHandler(
            filename = function() { paste0("volcano_plot_", Sys.Date(), ".png") },
            content = function(file) { ggsave(file, plot = volcano_plot(), width = 12, height = 10, dpi = 300) }
        )
    })
}






# volcanoUI <- function(id) {
#     
#     ns <- NS(id)
#     
#     sidebarLayout(
#         
#         sidebarPanel(
#             numericInput(ns("logfc_cutoff"), "Log2FC cutoff", 1),
#             selectInput(ns("p_metric"), "Significance metric:",
#                         choices = c("Adjusted p-value (padj)" = "padj",
#                                     "p-value" = "pvalue")),
#             numericInput(ns("p_cutoff"), "Cutoff", 0.05),
#             downloadButton(ns("download_plot"), "Download Plot")
#         ),
#         
#         mainPanel(
#             plotOutput(ns("volcano_plot"), height = "700px")
#         )
#         
#     )
#     
# }
# 
# 
# volcanoServer <- function(input, output, session, deseq_data) {
#     
#     volcano_plot <- reactive({
#         
#         df <- deseq_data()
#         
#         req(df)
#         req(all(c("log2FoldChange", "Geneid") %in% colnames(df)))
#         req(any(c("padj", "pvalue") %in% colnames(df)))  
#         
#         logfc_cutoff <- input$logfc_cutoff
#         p_cutoff <- input$p_cutoff
#         p_metric <- input$p_metric
#         
#         # Ensure chosen column exists
#         if (!(p_metric %in% colnames(df))) {
#             validate(need(FALSE, paste("Column", p_metric, "not found in data")))
#         }
#         
#         df <- df[!is.na(get(p_metric))]
#         df$y <- -log10(df[[p_metric]])
#         
#         # Annotate
#         df[, ann := fifelse(
#             get(p_metric) > p_cutoff, "Not significant",
#             fifelse(log2FoldChange > 0, "Up regulated", "Down regulated")
#         )]
#         
#         df$ann <- fifelse(
#             df[[p_metric]] <= p_cutoff & df$log2FoldChange > -1 & df$log2FoldChange < 1, 
#             paste0(df$ann, " (low)"),
#             df$ann
#         )
#         
#         df2 <- df[get(p_metric) <= p_cutoff & abs(log2FoldChange) >= logfc_cutoff]
#         df2 <- df2[order(abs(log2FoldChange), decreasing = TRUE)]
#         df2 <- df2[, head(.SD, 10), by = ann]
#         
#         ggplot(df, aes(x = log2FoldChange, y = y)) +
#             geom_point(aes(fill = ann), shape = 21, stroke = NA, size = 2, alpha = 0.5) +
#             geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linewidth = 0.3, linetype = "dashed") +
#             geom_hline(yintercept = -log10(p_cutoff), linewidth = 0.3, linetype = "dashed") +
#             geom_point(data = df2, aes(x = log2FoldChange, y = y, fill = ann), 
#                        shape = 21, stroke = 0.15, size = 2, color = "white") +
#             geom_text_repel(data = df2, aes(label = Geneid), max.overlaps = Inf,
#                             fontface = "bold", size = 2.5, bg.color = "white", bg.r = 0.05) +
#             scale_fill_manual(
#                 values = c(
#                     "Up regulated" = "#990000",
#                     "Up regulated (low)" = lighten("#990000", .5),
#                     "Down regulated" = "#004d99",
#                     "Down regulated (low)" = lighten("#004d99", .5),
#                     "Not significant" = "grey"
#                 ),
#                 breaks = c("Up regulated", "Not significant", "Down regulated"),
#                 guide = guide_legend(override.aes = list(size = 3, alpha = 1))
#             ) +
#             scale_x_continuous(breaks = c(-5, -2.5, -1, 0, 1, 2.5, 5),
#                                trans = scales::pseudo_log_trans()) +
#             scale_y_continuous(expand = c(0, 0),
#                                breaks = c(2, 5, 10, 20, 30, 40),
#                                trans = scales::pseudo_log_trans()) +
#             coord_cartesian(clip = "off") +
#             theme_minimal() +
#             theme(
#                 legend.title = element_blank(),
#                 legend.position = "bottom",
#                 axis.line = element_line(linewidth = .3, color = "black"),
#                 axis.ticks = element_line(linewidth = .3, color = "black"),
#                 panel.grid.minor = element_blank(),
#                 panel.grid.major = element_line(linewidth = .3, linetype = "dashed", color = "grey85"),
#                 plot.margin = margin(20, 20, 20, 20)
#             ) +
#             labs(x = "log2(Fold Change)", y = paste0("-log10(", p_metric, ")"))
#     })
#     
#     output$volcano_plot <- renderPlot({ volcano_plot() })
#     
#     output$download_plot <- downloadHandler(
#         
#         filename = function() {
#             paste0("volcano_plot_", Sys.Date(), ".png")
#         },
#         
#         content = function(file) {
#             ggsave(file, plot = volcano_plot(), width = 12, height = 10, dpi = 300)
#         }
#         
#     )
# }
# 
