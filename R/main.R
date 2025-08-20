
source("libraries.R")
source("upload_module.R")
source("volcano_module.R")
source("pca_module.R")   


ui <- navbarPage(
    
    title = "Seq2Viz",
    
    # theme = bs_theme(bootswatch = "flatly"),
    
    theme = bs_theme(
        preset = "cosmo",
        bg = "#F3F6FA",
        fg = "#386375",
        base_font = font_google("Jost"),
        "navbar-bg" = "#2b5769",   # dark blue/grey for navbar background
        "navbar-fg" = "#F3F6FA"    # white text/icons on the navbar
    ),

    # Welcome tab
    # tabPanel(
    #     "Home",
    #     fluidPage(
    #         titlePanel("Welcome"),
    #         fluidRow(
    #             column(
    #                 width = 8,
    #                 p("This app allows you to explore your data through interactive plots."),
    #                 p("You can upload your files on the first page, then navigate to:"),
    #                 tags$ul(
    #                     tags$li("Volcano plot tab to visualize differential expression."),
    #                     tags$li("PCA plot tab to explore sample clustering.")
    #                 ),
    #                 p("Start by uploading your data on the Upload tab!")
    #             )
    #         )
    #     )
    # ),
    
    # # Welcome tab
    # tabPanel(
    #     "Home",
    #     fluidPage(
    #         # titlePanel("Welcome"),
    #         fluidRow(
    #             column(
    #                 width = 8,
    #                 HTML("
    #     <div style='max-width: 800px; padding: 20px;'>
    #       <h3 style='color: #004164;'>Welcome</h3>
    #       <h5 style='color: #326286;'>
    #         This app allows you to explore RNA-seq data through interactive visualizations.
    #       </h5>
    # 
    #       <br>
    # 
    #       <h3 style='color: #004164;'>How to use this app</h3>
    #       <h5 style='color: #326286;'>
    #         You can upload your RNA-seq datasets on the Upload tab, then explore your data using the following interactive tools:
    #       </h5>
    #       <ul>
    #         <li>Volcano plot tab: visualize differential expression and highlight significant genes.</li>
    #         <li>PCA plot tab: explore sample clustering and relationships between samples.</li>
    #       </ul>
    #       <h5 style='color: #326286;'>
    #         Start by uploading your files on the Upload tab to get started!
    #       </h5>
    # 
    #       <br>
    # 
    #       <h3 style='color: #004164;'>Feedback & Contribution</h3>
    #       <h5 style='color: #326286;'>
    #         Your input helps improve the app. Report bugs or suggest new features via 
    #         <a href='https://github.com/your-repo/issues' style='color: #004164;'>GitHub Issues</a> 
    #         or contact the development team via 
    #         <a href='mailto:your-email@example.com' style='color: #004164;'>email</a>.
    #       </h5>
    # 
    #       <br>
    # 
    #       <h3 style='color: #004164;'>License</h3>
    #       <h5 style='color: #326286;'>
    #         This app is licensed under the MIT license. You are free to use and modify it, as long as the source is cited.
    #       </h5>
    #     </div>
    #     ")
    #             )
    #         )
    #     )
    # ),
    # 
    
    tabPanel(
        "Home",
        fluidPage(
            # titlePanel("Welcome to RNA-seq Explorer"),
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
                        h5("Begin by uploading your RNA-seq datasets on the Upload tab. Then navigate through the following features:"),
                        tags$ul(
                            tags$li("Volcano Plot: Visualize differential gene expression."),
                            tags$li("PCA Plot: Explore sample clustering and relationships.")
                        ),
                        style = "margin-bottom: 1rem; border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem;"
                    ),
                    
                    # Feedback Card
                    card(
                        title = "Feedback & Contribution",
                        h3("Feedback & Contribution", style = "color:#2b5769;"),
                        p(HTML("Report bugs or suggest new features via <a href='https://github.com/your-repo/issues'>GitHub Issues</a> or contact us via <a href='mailto:your-email@example.com'>email</a>.")),
                        style = "margin-bottom: 1rem; border-radius: 1rem; box-shadow: 2px 2px 10px rgba(0,0,0,0.1); padding: 1rem;"
                    ),
                    
                    # License Card
                    card(
                        title = "License",
                        h3("License", style = "color:#2b5769;"),
                        p("This app is licensed under the MIT license. You are free to use and modify it, as long as the source is cited."),
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
