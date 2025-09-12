
pcaUI <- function(id) {
    
    ns <- NS(id)
    
    sidebarLayout(
        
        sidebarPanel(
            
            uiOutput(ns("group")),
            uiOutput(ns("axis")), 
            downloadButton(ns("download_plot"), "Download Plot")
            
        ),
        
        mainPanel(
            plotOutput(ns("pca_plot"), height = "700px")
        )
        
    )
}

pcaServer <- function(input, output, session, meta_data, counts_data, deseq_data) {
    
    data_list <- reactive({
        
        req(meta_data(), counts_data(), deseq_data())
        meta <- meta_data()
        counts <- counts_data()
        deseq <- deseq_data()
        
        # Rename Geneid to gene_name if needed (match column name)
        if ("Geneid" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) {
            setnames(counts, "Geneid", "gene_name")
        }
        
        # --- filter significant genes ---
        sig_results <- deseq[which((padj <= 0.05) & abs(log2FoldChange) >= 2)]
        sig_genes   <- sig_results$Geneid
        counts <- counts[gene_name %in% sig_genes]

        list(meta = meta, counts = counts)
    })
    
    output$group <- renderUI({
        
        req(data_list())
        
        meta_cols <- names(data_list()$meta)
        
        selectInput(
            session$ns("group_col"), "Select grouping variable",
            choices = setdiff(meta_cols, "sampleID"),
            selected = if ("Group1" %in% meta_cols) "Group1" else meta_cols[1]
        )
    })
    
    # output$axis <- renderUI({
    #     
    #     req(data_list())
    #     
    #     meta_cols <- names(data_list()$meta)
    #     
    #     selectInput(
    #         session$ns("axis"), "Select axis variable",
    #         choices = 
    #     )
    # })
    
    output$axis <- renderUI({
        res <- pca_results()
        req(res)
        pcs <- names(res$pca_dt)
        pcs <- pcs[grepl("^PC\\d+$", pcs)]
        
        tagList(
            selectInput(session$ns("x_pc"), "X axis (PC)", choices = pcs,
                        selected = if ("PC1" %in% pcs) "PC1" else pcs[1]),
            selectInput(session$ns("y_pc"), "Y axis (PC)", choices = pcs,
                        selected = if ("PC2" %in% pcs) "PC2" else pcs[min(2, length(pcs))]),
            numericInput(session$ns("x_min"), "X min", value = -55),
            numericInput(session$ns("x_max"), "X max", value =  55),
            numericInput(session$ns("y_min"), "Y min", value = -55),
            numericInput(session$ns("y_max"), "Y max", value =  55)
        )
    })
    
    
    pca_results <- reactive({
        
        data <- data_list()
        meta <- data$meta
        counts <- data$counts
        group_col <- input$group_col
        
        
        validate(
            need("sampleID" %in% colnames(meta), "Sample metadata must include a 'sampleID' column."),
            need("gene_name" %in% colnames(counts), "Counts table must include a 'gene_name' column."),
            need(group_col %in% colnames(meta), "Selected grouping variable not found in metadata.")
        )
        
        # Prepare counts matrix: genes as rows, samples as columns
        gene_names <- counts$gene_name
        count_values <- as.matrix(counts[, setdiff(colnames(counts), "gene_name"), with = FALSE])
        rownames(count_values) <- gene_names
        
        # Transpose so samples are rows, genes are columns
        counts_mat <- t(count_values)
        
        # Match samples in metadata and counts
        common_samples <- intersect(rownames(counts_mat), meta$sampleID)
        counts_mat <- counts_mat[common_samples, , drop = FALSE]
        meta <- meta[match(common_samples, meta$sampleID)]
        
        # Run PCA
        pca <- prcomp(counts_mat, center = TRUE, scale. = TRUE)
        pca_dt <- as.data.table(pca$x, keep.rownames = "sampleID")
        
        # Merge PCA results with metadata
        pca_dt <- merge(pca_dt, meta, by = "sampleID", all.x = TRUE)
        
        list(pca_dt = pca_dt, pca = pca)
    })
    
    plot_pca <- reactive({
        res <- pca_results()
        pca_dt <- res$pca_dt
        pca <- res$pca
        group_col <- input$group_col
        
        # selected axes (fallbacks if UI not ready yet)
        xpc <- req(if (!is.null(input$x_pc)) input$x_pc else "PC1")
        ypc <- req(if (!is.null(input$y_pc)) input$y_pc else "PC2")
        
        # explained variance labels for selected PCs
        prop_var <- summary(pca)$importance[2, ]  # names like "PC1","PC2",...
        # make sure named; if not, set names
        if (is.null(names(prop_var))) names(prop_var) <- paste0("PC", seq_along(prop_var))
        x_lab <- paste0(xpc, " (", round(100 * (prop_var[[xpc]] %||% 0), 2), "%)")
        y_lab <- paste0(ypc, " (", round(100 * (prop_var[[ypc]] %||% 0), 2), "%)")
        
        label_col <- if ("patientID" %in% colnames(pca_dt)) "patientID" else "sampleID"
        
        # color palettes by group (unchanged logic)
        groups <- unique(pca_dt[[group_col]])
        groups <- groups[!is.na(groups)]
        n_groups <- length(groups)
        base_fill_colors <- c("#990000", "#004d99")
        base_color_colors <- c("#990000", "#004d99")
        pal_fill <- colorspace::lighten(
            grDevices::colorRampPalette(base_fill_colors)(max(1, n_groups)), amount = 0.25
        )
        pal_color <- colorspace::darken(
            grDevices::colorRampPalette(base_color_colors)(max(1, n_groups)), amount = 0.25
        )
        names(pal_fill) <- groups
        names(pal_color) <- groups
        
        # dynamic plot ranges (avoid fixed -55..55 which may crop other PCs)
        xr <- range(pca_dt[[xpc]], na.rm = TRUE)
        yr <- range(pca_dt[[ypc]], na.rm = TRUE)
        
        ggplot(pca_dt, aes_string(x = xpc, y = ypc, fill = group_col, color = group_col)) +
            ggforce::geom_mark_circle(alpha = 0.10, expand = unit(1.5, "mm")) +
            geom_point(shape = 21, size = 3, stroke = 0.25) +
            ggrepel::geom_text_repel(aes_string(label = label_col),
                                     fontface = "bold", size = 3,
                                     bg.color = "white", bg.r = 0.05) +
            scale_fill_manual(values = pal_fill, na.value = "grey") +
            scale_color_manual(values = pal_color, na.value = "grey") +
            scale_x_continuous(limits = c(input$x_min, input$x_max)) +
            scale_y_continuous(limits = c(input$y_min, input$y_max)) +
            theme_minimal() +
            theme(
                legend.position = "bottom",
                panel.grid = element_blank(),
                axis.line = element_line(lineend = "round"),
                axis.ticks = element_line(lineend = "round"),
                plot.margin = margin(20, 20, 20, 20)
            ) +
            labs(x = x_lab, y = y_lab)
    })
    
    # plot_pca <- reactive({
    #     
    #     res <- pca_results()
    #     pca_dt <- res$pca_dt
    #     pca <- res$pca
    #     group_col <- input$group_col
    #     axis <- input$axis
    #     
    #     prop_var <- summary(pca)$importance[2, ]
    #     pc1_lab <- paste0("PC1 (", round(prop_var[1] * 100, 2), "%)")
    #     pc2_lab <- paste0("PC2 (", round(prop_var[2] * 100, 2), "%)")
    #     
    #     label_col <- if ("patientID" %in% colnames(pca_dt)) "patientID" else "sampleID"
    #     
    #     # Get unique groups and remove NA
    #     groups <- unique(pca_dt[[group_col]])
    #     groups <- groups[!is.na(groups)]
    #     
    #     n_groups <- length(groups)
    #     base_fill_colors <- c("#990000", "#004d99")  
    #     base_color_colors <- c("#990000", "#004d99") 
    #     
    #     
    #     pal_fill <- colorspace::lighten(
    #         grDevices::colorRampPalette(base_fill_colors)(n_groups),
    #         amount = 0.25
    #     )
    #     pal_color <- colorspace::darken(
    #         grDevices::colorRampPalette(base_color_colors)(n_groups),
    #         amount = 0.25
    #     )
    #     
    #     names(pal_fill) <- groups
    #     names(pal_color) <- groups
    #     
    #     ggplot(pca_dt, aes_string(x = "PC1", y = "PC2", fill = group_col, color = group_col)) +
    #         
    #         ggforce::geom_mark_circle(alpha = 0.1, expand = unit(1.5, "mm")) +
    #         
    #         geom_point(shape = 21, size = 3, stroke = 0.25) +
    #         
    #         ggrepel::geom_text_repel(aes_string(label = label_col), fontface = "bold", size = 3, bg.color = "white", bg.r = 0.05) +
    #         
    #         scale_fill_manual(values = pal_fill, na.value = "grey") +
    #         
    #         scale_color_manual(values = pal_color, na.value = "grey") +
    #         
    #         scale_x_continuous(limits = c(-55, 55)) +
    #         scale_y_continuous(limits = c(-55, 55)) +
    #         
    #         theme_minimal() +
    #         
    #         theme(
    #             legend.position = "bottom",
    #             panel.grid = element_blank(),
    #             axis.line = element_line(lineend = "round"),
    #             axis.ticks = element_line(lineend = "round"),
    #             plot.margin = margin(20, 20, 20, 20)
    #         ) +
    #         
    #         labs(x = pc1_lab, y = pc2_lab)
    # })
    
    
    output$pca_plot <- renderPlot({
        plot_pca()
    })
    
    output$download_plot <- downloadHandler(
        filename = function() {
            paste0("PCA_plot_", Sys.Date(), ".png")
        },
        content = function(file) {
            ggsave(file, plot = plot_pca(), width = 10, height = 10, dpi = 300)
        }
    )
    
    
    
}
