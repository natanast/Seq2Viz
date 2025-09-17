

homeUI <- function(id) {
    ns <- NS(id)
    
    fluidPage(
        fluidRow(
            column(
                width = 12,
                
                # Welcome Card
                card(
                    title = "Welcome",
                    h4("Welcome", style = "color:#2b5769;"),
                    h6("Explore your RNA-seq data through interactive visualizations.", style = "color:#386375;"),
                    style = "margin-bottom: 1rem; border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem;"
                ),
                
                # How to Use Card 
                card(
                    title = "How to Use This App",
                    h4("How to Use This App", style = "color:#2b5769;"),
                    h5("Begin by uploading your RNA-seq datasets on the Upload tab (DESeq2 results, metadata, counts).", style = "color:#386375;"),
                    h5("File Requirements", style = "color:#2b5769; margin-top:0rem; margin-bottom:0rem;"),
                    tags$ul(
                        tags$li(HTML("<b>DESeq2 results file:</b> must contain <code>Geneid</code>, <code>pvalue</code>, and <code>padj</code> columns."), style = "color:#386375;"),
                        tags$li(HTML("<b>Gene counts file:</b> must contain a <code>gene_name</code> column."), style = "color:#386375;"),
                        tags$li(HTML("<b>Sample metadata file:</b> must contain <code>sampleID</code>, <code>patientID</code>, and <code>Group1</code> columns."), style = "color:#386375;")
                    ),
                    h5("Then navigate through the app features as follows:", style = "color:#386375; margin-top:0rem; margin-bottom:0rem;"),
                    tags$ol(
                        tags$li("PCA Plot: Explore sample clustering and relationships.", style = "color:#386375;"),
                        tags$li("Volcano Plot: Visualize differential gene expression.", style = "color:#386375;"),
                        tags$li("Heatmaps: Inspect gene expression patterns across samples.", style = "color:#386375;"),
                        # tags$li("GSEA analysis: Explore pathway enrichment.", style = "color:#386375;"),
                        # tags$li("ORA analysis: Run over-representation analysis.", style = "color:#386375;")
                    ),
                    h6("Tip: Hover over plots for detailed information. Adjust cutoffs and settings to filter results.", style = "color:#386375;"),
                    style = "margin-bottom: 1rem; border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem;"
                ),
                
                # Feedback Card
                card(
                    title = "Feedback & Contribution",
                    h4("Feedback & Contribution", style = "color:#2b5769;"),
                    h6(HTML("Report bugs or suggest new features via <a href='https://github.com/natanast/Seq2Viz/issues'>GitHub Issues</a>."),
                       style = "color:#386375;"),
                    style = "margin-bottom: 1rem; border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem;"
                )
                
            )
        )
    )
}

