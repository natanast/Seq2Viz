source("libraries.R")
source("upload.R")
source("de.R")
source("volcano.R")
source("pca.R")
source("home.R")
source("heatmap.R")



# Allow uploads up to 100 MB
options(shiny.maxRequestSize = 100*1024^2)


ui <- navbarPage(

    title = "Seq2Viz",

    theme = bs_theme(
        preset = "cosmo",
        bg = "#F3F6FA",
        fg = "#2b5769",
        base_font = font_google("Jost"),
        "navbar-fg" = "#F3F6FA"
    ),


    tabPanel(
        "Home",
        homeUI("home1")
    ),

    tabPanel(
        "Upload files",
         uploadUI("upload1")
    ),

    tabPanel(
        "Deferential exprassion",
        deUI("de")
    ),

    tabPanel(
        "PCA Plot",
         pcaUI("pca")
    ),

    tabPanel(
        "Volcano Plot",
         volcanoUI("volcano")
    ),

    tabPanel(
        "Heatmap",
         heatmapUI("heat1")
    )
)


server <- function(input, output, session) {
    
    # 1. Load Data
    data_list <- uploadServer("upload1")
    
    # 2. Run Analysis (Optional)
    de_results <- deserver("de", counts_data = data_list$counts, meta_data = data_list$metadata)
    
    # --- 3. THE SMART DATA SWITCH ---
    
    final_counts <- reactive({
        # Case A: User ran the DE tab analysis? -> Use calculated normalized counts
        if (!is.null(de_results()) && !is.null(de_results()$counts)) {
            return(de_results()$counts)
        } 
        # Case B: User uploaded a Normalized Counts file? -> Use that
        else if (!is.null(data_list$norm_counts())) {
            return(data_list$norm_counts())
        }
        # Case C: Fallback to raw counts (Visuals might be un-normalized)
        else {
            return(data_list$counts())
        }
    })
    
    final_deseq <- reactive({
        # Case A: User ran the DE tab analysis?
        if (!is.null(de_results()) && !is.null(de_results()$res)) {
            return(de_results()$res)
        }
        # Case B: User uploaded a DESeq results file?
        else {
            return(data_list$deseq())
        }
    })
    
    final_meta <- reactive({
        # Case A: User ran DE tab (factors might be updated)
        if (!is.null(de_results()) && !is.null(de_results()$meta)) {
            return(de_results()$meta)
        } 
        # Case B: Use uploaded metadata
        else {
            return(data_list$metadata())
        }
    })

    # 4. Run PCA with the chosen data
    pcaServer("pca", 
              meta_data = final_meta, 
              counts_data = final_counts, 
              deseq_data = final_deseq)
    
    
    # callModule(volcanoServer, "volcano", final_deseq)
    
    volcanoServer("volcano", deseq_data = final_deseq)


    # heatmapServer(
    #         "heat1",
    #         meta_data   = final_meta,
    #         counts_data = final_counts,
    #         deseq_data  = final_deseq
    #     )

    heatmapServer("heat1", meta_data = final_meta, counts_data = final_counts, deseq_data = final_deseq)
}

shinyApp(ui, server)




# server <- function(input, output, session) {
#     
#     data_list <- callModule(uploadServer, "upload1")
#     
#     # callModule(deserver, "de")
#     
#     # callModule(volcanoServer, "volcano", data_list$deseq)
#     # 
#     # callModule(pcaServer, "pca", data_list$metadata, data_list$counts, data_list$deseq)
#     # 
#     # heatmapServer(
#     #     "heat1",
#     #     meta_data   = data_list$metadata,
#     #     counts_data = data_list$counts,
#     #     deseq_data  = data_list$deseq
#     # )
#     # 
#     
# }
