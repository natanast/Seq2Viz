

heatmapUI <- function(id) {
    
    ns <- NS(id)
    
    tagList(
        
        sidebarLayout(
            
            sidebarPanel(
                
                h4("Data Parameters"),
                
                radioButtons(ns("mode"), "Mode:",
                             choices = c("Supervised (by Group)" = "supervised",
                                         "Unsupervised (cluster columns)" = "unsupervised"),
                             selected = "supervised"),
                
                radioButtons(ns("filter_column"), "Filter DE genes by:", 
                             choices = c("pvalue" = "pvalue", "padj" = "padj"), selected = "pvalue"),
                
                numericInput(ns("pval_thresh"), "p-value threshold", value = 0.05, min = 0, step = 0.01),
                numericInput(ns("lfc_thresh"), "abs(log2FoldChange) threshold (>)", value = 1, min = 0, step = 0.1),
                numericInput(ns("min_row_sum"), "Min gene total (rowSums)", value = 10, min = 0, step = 1),
                checkboxInput(ns("scale_rows"), "Scale rows (z-score)", value = TRUE),
                
                hr(),
                
                h4("Clustering Options"), # NEW SECTION
                selectInput(ns("dist_method"), "Distance Method",
                            choices = c("Pearson Correlation" = "pearson",
                                        "Euclidean" = "euclidean",
                                        "Spearman Correlation" = "spearman",
                                        "Manhattan" = "manhattan"),
                            selected = "pearson"),
                
                selectInput(ns("clust_method"), "Clustering Method",
                            choices = c("Ward.D2" = "ward.D2",
                                        "Complete" = "complete",
                                        "Average" = "average",
                                        "Single" = "single"),
                            selected = "ward.D2"),
                
                hr(),
                
                h4("Visual Parameters"),
                numericInput(ns("font_size"), "Font Size", value = 12, min = 6),
                checkboxInput(ns("show_row_names"), "Show Row Names", value = FALSE),
                
                conditionalPanel(condition = sprintf("input['%s'] == 'unsupervised'", ns("mode")),
                                 numericInput(ns("col_split"), "Column Split", value = 2, min = 1, step = 1)
                ),
                
                numericInput(ns("row_split"), "Row Split", value = 2, min = 1, step = 1),
                
                hr(),
                
                downloadButton(ns("download_heatmap"), "Download Heatmap (PNG)")
                
            ),
            
            mainPanel(
                plotOutput(ns("heatmap_plot"), height = "800px")
            )
            
        )
    )
}


heatmapServer <- function(id, meta_data, counts_data, deseq_data) {
    
    moduleServer(id, function(input, output, session) {
        
        # --- 1. Data Cleaning ---
        hm_data <- reactive({
            req(meta_data(), counts_data(), deseq_data())
            meta <- copy(meta_data())
            counts <- copy(counts_data())
            deseq <- copy(deseq_data())
            
            if ("Geneid" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) data.table::setnames(counts, "Geneid", "gene_name")
            if ("gene_id" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) data.table::setnames(counts, "gene_id", "gene_name")
            if (!("Sample" %in% colnames(meta)) && ("sampleID" %in% colnames(meta))) data.table::setnames(meta, "sampleID", "Sample")
            
            list(meta = meta, counts = counts, deseq = deseq)
        })
        
        # --- 2. Data Preparation ---
        heatmap_prep <- reactive({
            dd <- hm_data(); req(dd)
            meta <- dd$meta; counts <- dd$counts; deseq <- dd$deseq
            
            # Filter Significant Genes
            filter_col <- req(input$filter_column)
            d_sig <- deseq[get(filter_col) <= input$pval_thresh & abs(log2FoldChange) > input$lfc_thresh]
            
            gene_index <- sort(unique(d_sig$GeneID %||% d_sig$Geneid %||% d_sig$gene_name))
            if (length(gene_index) == 0) stop("No significant genes found.")
            
            # Filter Counts (Global)
            counts_filtered <- counts[gene_name %in% gene_index]
            if (nrow(counts_filtered) == 0) stop("No matching genes in counts.")
            
            # Prepare Matrix
            mm <- as.matrix(counts_filtered[, -1])
            rownames(mm) <- counts_filtered$gene_name
            
            # Row Sum Filter
            keep_rows <- which(rowSums(mm, na.rm=TRUE) >= input$min_row_sum)
            mm <- mm[keep_rows, , drop = FALSE]
            if (nrow(mm) == 0) stop("All genes filtered by row sum.")
            
            # Scale (Global Z-Score)
            if (isTRUE(input$scale_rows)) {
                mm <- t(scale(t(mm), center = TRUE, scale = TRUE))
                mm[is.na(mm)] <- 0
            }
            
            # Subset for View
            keep_samples <- intersect(meta$Sample, colnames(mm))
            if (length(keep_samples) == 0) stop("No overlapping samples found.")
            
            mm_view <- mm[, keep_samples]
            ann <- meta[Sample %in% keep_samples]
            ann <- ann[match(colnames(mm_view), ann$Sample)]
            
            if (input$mode == "supervised" && "Group" %in% colnames(ann)) {
                ann <- ann[order(ann$Group)]
                mm_view <- mm_view[, ann$Sample]
            }
            
            list(mat = mm_view, ann = ann)
        })
        
        # --- 3. Build Heatmap ---
        build_heatmap <- reactive({
            prep <- heatmap_prep(); req(prep)
            m <- prep$mat
            ann <- prep$ann
            
            set.seed(123)
            
            # Colors
            my_col <- colorRamp2(breaks = c(-4, -2, 0, 2, 4),
                                 colors = c('#00429d', '#73a2c6', 'grey96', '#f4777f', '#93003a'))
            
            # Dynamic Annotation
            ha <- NULL
            if ("Group" %in% colnames(ann)) {
                # Get groups in order (First level = Control/Reference)
                if (is.factor(ann$Group)) {
                    groups <- levels(ann$Group)
                } else {
                    groups <- unique(as.character(ann$Group))
                }
                
                # Define Ref (Blue) and Trt (Red)
                col_ref <- "#3a5cbc" 
                col_trt <- "#a62a17" 
                
                if (length(groups) == 2) {
                    final_colors <- c(col_ref, col_trt)
                } else {
                    final_colors <- colorRampPalette(c(col_ref, "#f4777f", col_trt))(length(groups))
                }
                
                pal <- setNames(final_colors, groups)
                
                ha <- HeatmapAnnotation(
                    Group = ann$Group, 
                    col = list(Group = pal),
                    simple_anno_size = unit(0.3, "cm"),
                    annotation_legend_param = list(title = "Group")
                )
            }
            
            # Splitting Logic
            if (input$mode == "supervised") {
                col_split_val <- if("Group" %in% colnames(ann)) ann$Group else NULL
                cluster_cols <- FALSE
            } else {
                col_split_val <- as.integer(input$col_split)
                cluster_cols <- TRUE
            }
            
            row_split_val <- as.integer(input$row_split %||% 1)
            if(row_split_val < 2) row_split_val <- NULL
            
            # Plot with Dynamic Clustering Parameters
            Heatmap(
                m,
                name = if (isTRUE(input$scale_rows)) "Z-score" else "Value",
                col = my_col,
                use_raster = TRUE,
                
                # --- NEW: Use inputs for clustering ---
                clustering_distance_rows = input$dist_method,
                clustering_distance_columns = input$dist_method,
                clustering_method_rows = input$clust_method,
                clustering_method_columns = input$clust_method,
                
                border = TRUE,
                rect_gp = gpar(col = "white", lwd = .25),
                column_split = col_split_val,
                cluster_columns = cluster_cols,
                row_split = row_split_val,
                top_annotation = ha,
                show_row_names = input$show_row_names,
                row_names_gp = gpar(fontsize = input$font_size * 0.8),
                column_names_gp = gpar(fontsize = input$font_size * 0.8),
                heatmap_legend_param = list(labels_gp = gpar(fontsize = input$font_size))
            )
        })
        
        # --- 4. Render & Download ---
        output$heatmap_plot <- renderPlot({
            
            ht <- build_heatmap(); req(ht)
            grid::grid.newpage()
            ComplexHeatmap::draw(ht, merge_legend = TRUE)
            
        }, res = 96)
        
        output$download_heatmap <- downloadHandler(
            
            filename = function() paste0("heatmap_", Sys.Date(), ".png"),
            
            content = function(file) {
                set.seed(123)
                ht <- build_heatmap()
                gr <- grid::grid.grabExpr(ComplexHeatmap::draw(ht, merge_legends = TRUE)) |> ggplotify::as.ggplot()
                ggsave(file, plot = gr, width = 10, height = 10, dpi = 600)
            }
            
        )
        
        observe({
            
            req(input$filter_column)
            updateNumericInput(session, "pval_thresh", label = if (input$filter_column == "pvalue") "p-value threshold" else "padj threshold")
            
        })
        
    })
}