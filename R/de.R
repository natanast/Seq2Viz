
deUI <- function(id) {
    
    ns <- NS(id)
    
    sidebarLayout(
        
        sidebarPanel(
            uiOutput(ns("group")),
            
            radioButtons(
                ns("filter_column"), "Filter DE genes by:",
                choices = c("padj" = "padj", "pvalue" = "pvalue"),
                selected = "padj"
            ),
            numericInput(ns("pval_thresh"), "padj threshold", value = 0.05, min = 0, step = 0.01),
            numericInput(ns("lfc_thresh"), "abs(log2FoldChange) threshold", value = 1, min = 0, step = 0.1),
            
            # label controls
            checkboxInput(ns("show_labels"), "Show sample labels", value = TRUE),
            uiOutput(ns("label_col")),
            uiOutput(ns("axis")),
            
            downloadButton(ns("download_plot"), "Download Plot")
        ),
        
        mainPanel(
            plotOutput(ns("pca_plot"), height = "700px")
        )
    )
}

deServer <- function(input, output, session, meta_data, counts_data, deseq_data) {
    
    ns <- session$ns
    
    # 1. Update numeric inputs labels
    observe({
        req(input$filter_column)
        new_label <- if (input$filter_column == "padj") "padj threshold" else "p-value threshold"
        updateNumericInput(session, "pval_thresh", label = new_label)
    })
    
    # 2. Data Preparation
    data_list <- reactive({
        
        req(meta_data(), counts_data(), deseq_data())
        meta <- meta_data()
        counts <- counts_data()
        deseq <- deseq_data()
        
        # Tolerant renaming for counts gene id col
        if ("Geneid" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) setnames(counts, "Geneid", "gene_name")
        if ("gene_id" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) setnames(counts, "gene_id", "gene_name")
        
        # Filtering Genes (Matches your script's logic)
        filter_col <- req(input$filter_column)
        pth <- req(input$pval_thresh)
        lfc_cut <- req(input$lfc_thresh)
        
        if (!(filter_col %in% colnames(deseq))) stop(sprintf("DE table does not contain column '%s'.", filter_col))
        if (!("log2FoldChange" %in% colnames(deseq))) stop("DE table must contain 'log2FoldChange' column.")
        
        d_sig <- deseq[ which( get(filter_col) <= pth & abs(log2FoldChange) >= lfc_cut ) ]
        gene_ids <- unique(d_sig$GeneID %||% d_sig$Geneid %||% d_sig$gene_id)
        
        # Filter Counts by GENE only (Keep ALL samples for now, just like the script)
        counts_filtered <- if (length(gene_ids) == 0) counts[0, ] else counts[ gene_name %in% gene_ids ]
        
        list(meta = meta, counts = counts_filtered)
        
    })
    
    # 3. UI Updates
    output$group <- renderUI({
        
        req(data_list())
        meta_cols <- names(data_list()$meta)
        selectInput(
            ns("group_col"), "Select grouping variable",
            choices = setdiff(meta_cols, "sampleID"),
            selected = if ("Group" %in% meta_cols) "Group" else setdiff(meta_cols, "sampleID")[1]
        )
        
    })
    
    output$label_col <- renderUI({
        
        req(data_list())
        meta_cols <- names(data_list()$meta)
        default <- if ("patientID" %in% meta_cols) "patientID" else if ("sampleID" %in% meta_cols) "sampleID" else meta_cols[1]
        selectInput(ns("label_col"), "Label column", choices = meta_cols, selected = default)
        
    })
    
    output$axis <- renderUI({
        
        res_ok <- tryCatch({ pca_results(); TRUE }, error = function(e) FALSE)
        pcs <- NULL
        if (res_ok) {
            res <- pca_results()
            req(res)
            # Only show PCs that actually exist in the result
            pcs <- colnames(res$pca$x)
        }
        if (is.null(pcs) || length(pcs) == 0) pcs <- c("PC1", "PC2", "PC3", "PC4")
        
        tagList(
            selectInput(ns("x_pc"), "X axis (PC)", choices = pcs, selected = if ("PC1" %in% pcs) "PC1" else pcs[1]),
            selectInput(ns("y_pc"), "Y axis (PC)", choices = pcs, selected = if ("PC2" %in% pcs) "PC2" else pcs[min(2, length(pcs))]),
            numericInput(ns("x_min"), "X min", value = -50),
            numericInput(ns("x_max"), "X max", value =  50),
            numericInput(ns("y_min"), "Y min", value = -50),
            numericInput(ns("y_max"), "Y max", value =  50)
        )
        
    })
    
    # 4. PCA Calculation 
    pca_results <- reactive({
        
        data <- data_list()
        meta <- data$meta
        counts <- data$counts
        
        validate(need(nrow(counts) > 0, "No genes passed the selected significance filters. Relax thresholds."))
        
        # --- Prepare Matrix for PCA ---
        gene_names <- counts$gene_name
        sample_cols <- setdiff(colnames(counts), "gene_name")
        
        # Create matrix of ALL samples in the counts file
        count_values <- as.matrix(counts[, sample_cols, with = FALSE])
        rownames(count_values) <- gene_names
        counts_mat <- t(count_values) # samples as rows, genes as columns
        
        # --- Run PCA on ALL samples ---
        pca <- prcomp(counts_mat, center = TRUE, scale. = TRUE)
        
        # Extract PCA coordinates
        pca_dt <- as.data.table(pca$x, keep.rownames = "sampleID")
        
        # Check intersection
        common_samples <- intersect(pca_dt$sampleID, meta$sampleID)
        validate(need(length(common_samples) > 0, "No matching sample names between counts and metadata."))
        
        # Filter the RESULT of the PCA, not the INPUT
        pca_dt <- pca_dt[sampleID %in% common_samples]
        
        # Merge with metadata
        pca_dt <- merge(pca_dt, meta, by = "sampleID", all.x = TRUE)
        
        list(pca_dt = pca_dt, pca = pca)
        
    })
    
    # 5. Plotting
    plot_pca <- reactive({
        
        res <- pca_results()
        pca_dt <- res$pca_dt
        pca <- res$pca
        
       
        req(input$group_col)
        group_col <- input$group_col
        
        # Safety check: ensure the column is actually in the data
        validate(need(group_col %in% colnames(pca_dt), 
                      paste0("Column '", group_col, "' not found in metadata.")))
        
        
        xpc <- req(if (!is.null(input$x_pc)) input$x_pc else "PC1")
        ypc <- req(if (!is.null(input$y_pc)) input$y_pc else "PC2")
        
        # Calculate variance explained
        prop_var <- summary(pca)$importance[2, ]
        if (is.null(names(prop_var))) names(prop_var) <- paste0("PC", seq_along(prop_var))
        
        idx_x <- as.integer(sub("^PC", "", xpc))
        idx_y <- as.integer(sub("^PC", "", ypc))
        
        pc_x_pct <- if (!is.na(idx_x) && idx_x <= length(prop_var)) prop_var[idx_x] else NA
        pc_y_pct <- if (!is.na(idx_y) && idx_y <= length(prop_var)) prop_var[idx_y] else NA
        
        x_lab <- if (!is.na(pc_x_pct)) paste0(xpc, " (", round(pc_x_pct * 100, 2), "%)") else xpc
        y_lab <- if (!is.na(pc_y_pct)) paste0(ypc, " (", round(pc_y_pct * 100, 2), "%)") else ypc
        
        label_col <- if (!is.null(input$label_col)) input$label_col else (if ("patientID" %in% colnames(pca_dt)) "patientID" else "sampleID")
        
        groups <- unique(pca_dt[[group_col]])
        groups <- groups[!is.na(groups)]
        n_groups <- length(groups)
        
        base_fill_colors <- c("#990000", "#004d99")
        base_color_colors <- c("#990000", "#004d99")
        
        pal_fill <- colorspace::lighten(grDevices::colorRampPalette(base_fill_colors)(max(1, n_groups)), amount = 0.25)
        pal_color <- colorspace::darken(grDevices::colorRampPalette(base_color_colors)(max(1, n_groups)), amount = 0.25)
        
        names(pal_fill) <- groups
        names(pal_color) <- groups
        
        x_limits <- c(input$x_min, input$x_max)
        y_limits <- c(input$y_min, input$y_max)
        
        p <- ggplot(pca_dt, aes_string(x = xpc, y = ypc, fill = group_col, color = group_col)) +
            ggforce::geom_mark_circle(alpha = 0.10, expand = unit(1.5, "mm")) +
            geom_point(shape = 21, size = 3, stroke = 0.25) +
            scale_fill_manual(values = pal_fill, na.value = "grey") +
            scale_color_manual(values = pal_color, na.value = "grey") +
            scale_x_continuous(limits = x_limits) +
            scale_y_continuous(limits = y_limits) +
            theme_minimal() +
            theme(
                legend.position = "bottom",
                panel.grid = element_blank(),
                axis.line = element_line(lineend = "round"),
                axis.ticks = element_line(lineend = "round"),
                plot.margin = margin(20, 20, 20, 20)
            ) +
            labs(x = x_lab, y = y_lab)
        
        if ( isTRUE(input$show_labels) ) {
            
            if (!(label_col %in% colnames(pca_dt))) label_col <- "sampleID"
            p <- p + geom_text_repel(
                aes_string(label = label_col),
                fontface = "bold", size = 3,
                bg.color = "white", bg.r = 0.05
            )
        }
        
        p
    })
    
    
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

