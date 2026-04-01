
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
    de_state <- deserver("de", counts_data = data_list$counts, meta_data = data_list$metadata)
    
    # --- 3. THE SMART DATA SWITCH ---

    active_subset <- reactive({
        de_subset <- de_state$active_subset()
        if (!is.null(de_subset)) {
            return(de_subset)
        }
        data_list$external_subset()
    })
    
    filtered_counts <- reactive({
        subset_info <- active_subset()
        if (!is.null(subset_info) && !is.null(subset_info$counts)) {
            return(subset_info$counts)
        }
        data_list$counts()
    })

    filtered_meta <- reactive({
        subset_info <- active_subset()
        if (!is.null(subset_info) && !is.null(subset_info$meta)) {
            return(subset_info$meta)
        }
        data_list$metadata()
    })

    subset_counts_to_meta <- function(counts_dt, meta_dt) {
        if (is.null(counts_dt) || is.null(meta_dt)) return(counts_dt)

        counts_dt <- copy(counts_dt)
        meta_dt <- copy(meta_dt)

        gene_candidates <- c("gene_name", "Geneid", "gene_id", "GeneID", "id", "ID", "gene", "Gene")
        gene_col <- intersect(gene_candidates, colnames(counts_dt))[1]
        if (is.na(gene_col)) gene_col <- colnames(counts_dt)[1]

        sample_col <- grep("sample|id", colnames(meta_dt), ignore.case = TRUE, value = TRUE)[1]
        if (is.na(sample_col)) sample_col <- colnames(meta_dt)[1]

        keep_cols <- c(gene_col, intersect(colnames(counts_dt), meta_dt[[sample_col]]))

        out <- counts_dt[, ..keep_cols]
        if ("gene_name" %in% colnames(out)) {
            setcolorder(out, c("gene_name", setdiff(colnames(out), "gene_name")))
        } else if (gene_col %in% colnames(out) && gene_col != colnames(out)[1]) {
            setcolorder(out, c(gene_col, setdiff(colnames(out), gene_col)))
        }
        out
    }

    pca_counts <- reactive({
        # Case A: User ran the DE tab analysis? -> Use calculated normalized counts
        if (!is.null(de_state$results()) && !is.null(de_state$results()$counts)) {
            return(de_state$results()$counts)
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

    heatmap_counts <- reactive({
        # Case A: User ran the DE tab analysis? -> Use calculated normalized counts
        if (!is.null(de_state$results()) && !is.null(de_state$results()$counts)) {
            return(de_state$results()$counts)
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
        if (!is.null(de_state$results()) && !is.null(de_state$results()$res)) {
            return(de_state$results()$res)
        }
        # Case B: User uploaded a DESeq results file?
        else {
            return(data_list$deseq())
        }
    })
    
    final_meta <- reactive({
        filtered_meta()
    })

    heatmap_display_meta <- reactive({
        meta <- filtered_meta()
        subset_info <- active_subset()

        if (is.null(meta) || is.null(subset_info)) {
            return(meta)
        }

        main_factor <- subset_info$main_factor
        if (!main_factor %in% colnames(meta)) {
            return(meta)
        }

        ordered_levels <- c(
            subset_info$target_level,
            subset_info$ref_level,
            setdiff(unique(as.character(meta[[main_factor]])), c(subset_info$target_level, subset_info$ref_level))
        )

        meta[[main_factor]] <- factor(as.character(meta[[main_factor]]), levels = ordered_levels)
        meta
    })

    full_meta <- reactive({
        if (!is.null(de_state$results()) && !is.null(de_state$results()$meta)) {
            return(de_state$results()$meta)
        }
        data_list$metadata()
    })

    # 4. Run PCA with the chosen data
    pcaServer("pca", 
              meta_data = final_meta, 
              counts_data = pca_counts, 
              deseq_data = final_deseq)
    
    
    # callModule(volcanoServer, "volcano", final_deseq)
    
    volcanoServer("volcano", deseq_data = final_deseq)


    # heatmapServer(
    #         "heat1",
    #         meta_data   = final_meta,
    #         counts_data = final_counts,
    #         deseq_data  = final_deseq
    #     )

    heatmapServer(
        "heat1",
        meta_data = full_meta,
        display_meta_data = heatmap_display_meta,
        counts_data = heatmap_counts,
        deseq_data = final_deseq
    )
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
