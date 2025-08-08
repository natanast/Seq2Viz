pca_ui <- function(id) {
    ns <- NS(id)
    tagList(
        sidebarLayout(
            sidebarPanel(
                fileInput(ns("metadata_file"), "Upload Sample Metadata (.xlsx, .csv, .tsv)"),
                fileInput(ns("counts_file"), "Upload Normalized Counts (.txt, .csv, .tsv)"),
                uiOutput(ns("group_selector"))
            ),
            mainPanel(
                plotOutput(ns("pca_plot"), height = "700px")
            )
        )
    )
}

pca_server <- function(id) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns
        
        data_list <- reactive({
            req(input$metadata_file, input$counts_file)
            
            # Load metadata
            meta_ext <- tools::file_ext(input$metadata_file$name)
            meta <- switch(
                meta_ext,
                xlsx = as.data.table(readxl::read_xlsx(input$metadata_file$datapath)),
                xls = as.data.table(readxl::read_xlsx(input$metadata_file$datapath)),
                csv = fread(input$metadata_file$datapath),
                tsv = fread(input$metadata_file$datapath, sep = "\t"),
                txt = fread(input$metadata_file$datapath, sep = "\t"),
                stop("Unsupported metadata file format")
            )
            
            # Load counts
            counts_ext <- tools::file_ext(input$counts_file$name)
            counts <- switch(
                counts_ext,
                csv = fread(input$counts_file$datapath),
                txt = fread(input$counts_file$datapath),
                tsv = fread(input$counts_file$datapath, sep = "\t"),
                stop("Unsupported counts file format")
            )
            
            # Rename Geneid to gene_name to match code expectation
            if ("Geneid" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) {
                setnames(counts, "Geneid", "gene_name")
            }
            
            list(meta = meta, counts = counts)
        })
        
        # Dynamically select grouping variable from metadata columns (except sampleID)
        output$group_selector <- renderUI({
            req(data_list())
            meta_cols <- names(data_list()$meta)
            selectInput(ns("group_col"), "Select grouping variable",
                        choices = setdiff(meta_cols, "sampleID"),
                        selected = "Group1")  # default to Group1 since you have Group1 and Group2
        })
        
        pca_results <- reactive({
            data <- data_list()
            meta <- data$meta
            counts <- data$counts
            
            validate(
                need("sampleID" %in% colnames(meta), "Sample metadata must include a 'sampleID' column."),
                need("gene_name" %in% colnames(counts), "Counts table must include a 'gene_name' column."),
                need(input$group_col %in% colnames(meta), "Selected grouping variable not found in metadata.")
            )
            
            # Transform counts: transpose genes to columns, samples to rows
            counts_t <- transpose(counts, make.names = "gene_name", keep.names = "sampleID")
            counts_mat <- counts_t[, -1, with = FALSE]
            rownames(counts_mat) <- counts_t$sampleID
            
            # PCA
            pca <- prcomp(counts_mat, scale. = TRUE, center = TRUE)
            pca_dt <- as.data.table(pca$x, keep.rownames = "sampleID")
            
            # Merge PCA results with metadata
            pca_dt <- merge(pca_dt, meta, by = "sampleID", all.x = TRUE)
            
            list(pca_dt = pca_dt, pca = pca)
        })
        
        output$pca_plot <- renderPlot({
            res <- pca_results()
            pca_dt <- res$pca_dt
            pca <- res$pca
            
            group_col <- input$group_col
            
            prop_var <- summary(pca)$importance[2, ]
            pc1_lab <- paste0("PC1 (", round(prop_var[1] * 100, 2), "%)")
            pc2_lab <- paste0("PC2 (", round(prop_var[2] * 100, 2), "%)")
            
            ggplot(pca_dt, aes_string(x = "PC1", y = "PC2", fill = group_col, color = group_col)) +
                ggforce::geom_mark_circle(alpha = 0.1, expand = unit(1.5, "mm")) +
                geom_point(shape = 21, size = 3, stroke = 0.25) +
                ggrepel::geom_text_repel(aes(label = patientID), fontface = "bold", size = 2, bg.color = "white", bg.r = 0.05) +
                scale_fill_manual(values = c("BIRC3_c1639del" = "#990000", "empty_backbone" = "#004d99") |> colorspace::lighten(0.25)) +
                scale_color_manual(values = c("BIRC3_c1639del" = "#990000", "empty_backbone" = "#004d99") |> colorspace::darken(0.25)) +
                theme_minimal() +
                theme(legend.position = "bottom",
                      panel.grid = element_blank(),
                      axis.line = element_line(lineend = "round"),
                      axis.ticks = element_line(lineend = "round"),
                      plot.margin = margin(20, 20, 20, 20)) +
                labs(x = pc1_lab, y = pc2_lab)
        })
    })
}
