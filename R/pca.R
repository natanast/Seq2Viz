
pcaUI <- function(id) {
    ns <- NS(id)
    
    sidebarLayout(
        sidebarPanel(
            
            uiOutput(ns("group")),
            
            # filtering controls for DE genes
            radioButtons(
                ns("filter_column"), "Filter DE genes by:",
                choices = c("padj" = "padj", "pvalue" = "pvalue"),
                selected = "padj"
            ),
            
            # threshold for chosen filter 
            numericInput(ns("pval_thresh"), "padj threshold", value = 0.05, min = 0, step = 0.01),
            
            # absolute log2FC cutoff (user-changeable)
            numericInput(ns("lfc_thresh"), "abs(log2FoldChange) threshold", value = 1, min = 0, step = 0.1),
            
            uiOutput(ns("axis")),
            
            downloadButton(ns("download_plot"), "Download Plot")
        ),
        
        mainPanel(
            plotOutput(ns("pca_plot"), height = "700px")
        )
    )
}

pcaServer <- function(input, output, session, meta_data, counts_data, deseq_data) {
    
    ns <- session$ns
    
    # observe filter_column and update numeric input label accordingly
    observe({
        req(input$filter_column)
        new_label <- if (input$filter_column == "padj") "padj threshold" else "p-value threshold"
        updateNumericInput(session, "pval_thresh", label = new_label)
    })
    
    data_list <- reactive({
        
        req(meta_data(), counts_data(), deseq_data())
        meta <- meta_data()
        counts <- counts_data()
        deseq <- deseq_data()
        
        # tolerant renaming for counts gene id col
        if ("Geneid" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) {
            setnames(counts, "Geneid", "gene_name")
        }
        if ("gene_id" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) {
            setnames(counts, "gene_id", "gene_name")
        }
        
        # choose filter column and thresholds
        filter_col <- req(input$filter_column)    # "padj" or "pvalue"
        pth <- req(input$pval_thresh)
        lfc_cut <- req(input$lfc_thresh)
        
        # check deseq has required columns
        if (!(filter_col %in% colnames(deseq))) {
            stop(sprintf("DE table does not contain column '%s'.", filter_col))
        }
        if (!("log2FoldChange" %in% colnames(deseq))) {
            stop("DE table must contain 'log2FoldChange' column.")
        }
        
        # filter significant genes using the selected column and lfc cutoff
        # use data.table syntax safely
        d_sig <- deseq[ which( get(filter_col) <= pth & abs(log2FoldChange) >= lfc_cut ) ]
        
        # tolerant gene id extraction
        gene_ids <- unique(d_sig$GeneID %||% d_sig$Geneid %||% d_sig$gene_id)
        if (length(gene_ids) == 0) {
            # no genes passed filter -> keep counts empty (user will see message later)
            counts_filtered <- counts[0, ]
        } else {
            counts_filtered <- counts[ gene_name %in% gene_ids ]
        }
        
        list(meta = meta, counts = counts_filtered)
    })
    
    
    output$group <- renderUI({
        
        req(data_list())
        
        meta_cols <- names(data_list()$meta)
        
        selectInput(
            ns("group_col"), "Select grouping variable",
            choices = setdiff(meta_cols, "sampleID"),
            selected = if ("Group1" %in% meta_cols) "Group1" else setdiff(meta_cols, "sampleID")[1]
        )
    })
    
    
    output$axis <- renderUI({
        
        # need pca results to know available PCs — but if none yet, show defaults
        res_ok <- tryCatch({ pca_results(); TRUE }, error = function(e) FALSE)
        
        pcs <- NULL
        if (res_ok) {
            res <- pca_results()
            req(res)
            pcs <- names(res$pca_dt)
            pcs <- pcs[grepl("^PC\\d+$", pcs)]
        }
        
        if (is.null(pcs) || length(pcs) == 0) {
            # fallback choices
            pcs <- c("PC1", "PC2", "PC3", "PC4")
        }
        
        tagList(
            selectInput(ns("x_pc"), "X axis (PC)", choices = pcs,
                        selected = if ("PC1" %in% pcs) "PC1" else pcs[1]),
            
            selectInput(ns("y_pc"), "Y axis (PC)", choices = pcs,
                        selected = if ("PC2" %in% pcs) "PC2" else pcs[min(2, length(pcs))]),
            
            numericInput(ns("x_min"), "X min", value = -50),
            numericInput(ns("x_max"), "X max", value =  50),
            numericInput(ns("y_min"), "Y min", value = -50),
            numericInput(ns("y_max"), "Y max", value =  50)
        )
    })
    
    
    pca_results <- reactive({
        
        data <- data_list()
        meta <- data$meta
        counts <- data$counts
        group_col <- input$group_col
        
        # validations
        validate(
            need("sampleID" %in% colnames(meta), "Sample metadata must include a 'sampleID' column."),
            need("gene_name" %in% colnames(counts), "Counts table must include a 'gene_name' column."),
            need(group_col %in% colnames(meta), "Selected grouping variable not found in metadata."),
            need(nrow(counts) > 0, "No genes passed the selected significance filters. Relax thresholds.")
        )
        
        # Prepare counts matrix: genes as rows, samples as columns
        gene_names <- counts$gene_name
        # ensure we only keep actual sample columns (not accidental non-sample cols)
        sample_cols <- setdiff(colnames(counts), "gene_name")
        if (length(sample_cols) == 0) stop("Counts table contains no sample columns.")
        
        count_values <- as.matrix(counts[, sample_cols, with = FALSE])
        rownames(count_values) <- gene_names
        
        # Transpose so samples are rows, genes are columns
        counts_mat <- t(count_values)
        
        # Match samples in metadata and counts
        common_samples <- intersect(rownames(counts_mat), meta$sampleID)
        validate(need(length(common_samples) > 0, "No matching sample names between counts and metadata."))
        counts_mat <- counts_mat[common_samples, , drop = FALSE]
        meta <- meta[match(common_samples, meta$sampleID), , drop = FALSE]
        
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
        
        # selected axes
        xpc <- req(if (!is.null(input$x_pc)) input$x_pc else "PC1")
        ypc <- req(if (!is.null(input$y_pc)) input$y_pc else "PC2")
        
        # proportion of variance for selected PCs
        prop_var <- summary(pca)$importance[2, ]
        if (is.null(names(prop_var))) names(prop_var) <- paste0("PC", seq_along(prop_var))
        
        idx_x <- as.integer(sub("^PC", "", xpc))
        idx_y <- as.integer(sub("^PC", "", ypc))
        
        pc_x_pct <- if (!is.na(idx_x) && idx_x <= length(prop_var)) prop_var[idx_x] else NA
        pc_y_pct <- if (!is.na(idx_y) && idx_y <= length(prop_var)) prop_var[idx_y] else NA
        
        x_lab <- if (!is.na(pc_x_pct)) paste0(xpc, " (", round(pc_x_pct * 100, 2), "%)") else xpc
        y_lab <- if (!is.na(pc_y_pct)) paste0(ypc, " (", round(pc_y_pct * 100, 2), "%)") else ypc
        
        label_col <- if ("patientID" %in% colnames(pca_dt)) "patientID" else "sampleID"
        
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
        
        # dynamic plot ranges: use user inputs
        x_limits <- c(input$x_min, input$x_max)
        y_limits <- c(input$y_min, input$y_max)
        
        ggplot(pca_dt, aes_string(x = xpc, y = ypc, fill = group_col, color = group_col)) +
            
            geom_mark_circle(alpha = 0.10, expand = unit(1.5, "mm")) +
            
            geom_point(shape = 21, size = 3, stroke = 0.25) +
            
            geom_text_repel(aes_string(label = label_col),
                                     fontface = "bold", size = 3,
                                     bg.color = "white", bg.r = 0.05) +
            
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

