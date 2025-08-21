
source("libraries.R")
source("upload_module.R")
source("volcano_module.R")
source("pca_module.R")   


ui <- navbarPage(
    
    title = "Seq2Viz",
    
    theme = bs_theme(
        preset = "cosmo",
        bg = "#F3F6FA",
        fg = "#2b5769",
        base_font = font_google("Jost"),
        # "navbar-bg" = "#2b5769",   
        "navbar-fg" = "#F3F6FA"
    ),


    tabPanel(
        "Home",
        fluidPage(
            
            fluidRow(
                column(
                    width = 12,
                    
                    # Welcome Card
                    card(
                        title = "Welcome",
                        h3("Welcome", style = "color:#2b5769;"),
                        h5("Explore your RNA-seq data through interactive visualizations.", style = "color:#386375;"),
                        style = "margin-bottom: 1rem; border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem;"
                    ),
                    
                    # How to Use Card
                    card(
                        title = "How to Use This App",
                        h3("How to Use This App", style = "color:#2b5769;"),
                        h5("Begin by uploading your RNA-seq datasets on the Upload tab.",
                           style = "color:#386375;"),
                        h5("Then navigate through the following features:", style = "color:#386375;"),
                        tags$ul(
                            tags$li("PCA Plot: Explore sample clustering and relationships.", style = "color:#386375;"),
                            tags$li("Volcano Plot: Visualize differential gene expression.", style = "color:#386375;"),
                            tags$li("Heatmaps", style = "color:#386375;"),
                            tags$li("GSEA analysis", style = "color:#386375;"),
                            tags$li("ORA analysis", style = "color:#386375;"),
                        ),
                        style = "margin-bottom: 1rem; border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem;"
                    ),
                    
                    # Feedback Card
                    card(
                        title = "Feedback & Contribution",
                        h3("Feedback & Contribution", style = "color:#2b5769;"),
                        h5(HTML("Report bugs or suggest new features via <a href='https://github.com/your-repo/issues'>GitHub Issues</a> or contact us via <a href='mailto:your-email@example.com'>email</a>."),
                           style = "color:#386375;"),
                        style = "margin-bottom: 1rem; border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem;"
                    ),
                    
                    # License Card
                    card(
                        title = "License",
                        h3("License", style = "color:#2b5769;"),
                        h5("This app is licensed under the MIT license. You are free to use and modify it, as long as the source is cited.", 
                           style = "color:#386375;"),
                        style = "margin-bottom: 1rem; border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem;"
                    )
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
    callModule(pcaServer, "pca1", data_list$metadata, data_list$counts, data_list$deseq)
    
}


shinyApp(ui, server)
