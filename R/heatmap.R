

heatmapUI <- function(id) {
    ns <- NS(id)
    tagList(
        sidebarLayout(
            sidebarPanel(
                radioButtons(ns("mode"), "Mode:",
                             choices = c("Supervised (by Group)" = "supervised",
                                         "Unsupervised (cluster columns)" = "unsupervised"),
                             selected = "supervised"),
                radioButtons(
                    ns("filter_column"), 
                    "Filter DE genes by:", 
                    choices = c("pvalue" = "pvalue", "padj" = "padj"), 
                    selected = "pvalue"
                ),
                numericInput(ns("pval_thresh"), "p-value threshold", value = 0.05, min = 0, step = 0.01),
                numericInput(ns("lfc_thresh"), "abs(log2FoldChange) threshold", value = 1, min = 0, step = 0.1),
                numericInput(ns("min_row_sum"), "Min gene total across samples (rowSums)", value = 10, min = 0, step = 1),
                checkboxInput(ns("scale_rows"), "Scale rows (z-score)", value = TRUE),
                conditionalPanel(condition = sprintf("input['%s'] == 'unsupervised'", ns("mode")),
                                 numericInput(ns("col_split"), "column_split (numeric)", value = 2, min = 1, step = 1)
                ),
                numericInput(ns("row_split"), "Row split", value = 2, min = 1, step = 1),
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
        
        hm_data <- reactive({
            req(meta_data(), counts_data(), deseq_data())
            meta <- copy(meta_data())
            counts <- copy(counts_data())
            deseq <- copy(deseq_data())
            
            # tolerant renaming
            if ("Geneid" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) {
                data.table::setnames(counts, "Geneid", "gene_name")
            }
            if ("gene_id" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) {
                data.table::setnames(counts, "gene_id", "gene_name")
            }
            if (!("Sample" %in% colnames(meta)) && ("sampleID" %in% colnames(meta))) {
                data.table::setnames(meta, "sampleID", "Sample")
            }
            
            list(meta = meta, counts = counts, deseq = deseq)
        })
        
        heatmap_prep <- reactive({
            dd <- hm_data(); req(dd)
            meta <- dd$meta; counts <- dd$counts; deseq <- dd$deseq
            
            # validate
            req("Sample" %in% colnames(meta), "metadata must contain a 'Sample' column")
            req("gene_name" %in% colnames(counts), "counts must contain 'gene_name' column")
            req(("pvalue" %in% colnames(deseq) && "log2FoldChange" %in% colnames(deseq)),
                "deseq must contain 'pvalue' and 'log2FoldChange' columns")
            
            # Determine which column to filter on (pvalue or padj)
            filter_col <- req(input$filter_column)  # user choice: "pvalue" or "padj"
            if (!(filter_col %in% colnames(deseq))) stop(paste0("Column '", filter_col, "' not found in DESeq2 results."))
            
            # Filter DE genes
            pth <- req(input$pval_thresh)
            lfc <- req(input$lfc_thresh)
            d_sig <- deseq[get(filter_col) <= pth & abs(log2FoldChange) >= lfc]
            
            gene_index <- sort(unique(d_sig$GeneID %||% d_sig$Geneid %||% d_sig$gene_id))
            if (length(gene_index) == 0) stop("No significant genes found with the chosen thresholds.")
            
            counts_sub <- counts[gene_name %in% gene_index]
            if (nrow(counts_sub) == 0) stop("No matching genes found in counts for DE gene list.")
            
            # Select metadata rows corresponding to columns in counts_sub
            sample_cols <- setdiff(colnames(counts_sub), "gene_name")
            ann <- meta[Sample %in% sample_cols]
            
            # For supervised, ensure Group exists and set factor levels like your script
            if (input$mode == "supervised") {
                if (!("Group" %in% colnames(ann))) stop("Supervised mode requires 'Group' column in metadata.")
                # enforce factor levels if both known
                if (all(c("control", "ibrutinib") %in% unique(as.character(ann$Group)))) {
                    ann$Group <- factor(ann$Group, levels = c("control", "ibrutinib"))
                } else {
                    ann$Group <- as.factor(ann$Group)
                }
                # order samples by Group
                ann <- ann[order(ann$Group)]
            } else {
                # Unsupervised: keep metadata ordering as-is (or by Sample)
                ann <- ann[order(ann$Sample)]
            }
            
            # Reorder counts columns to match ann$Sample
            sel_samples <- intersect(ann$Sample, sample_cols)
            if (length(sel_samples) == 0) stop("No sample names overlap between metadata$Sample and counts columns.")
            mm <- as.matrix(counts_sub[, sel_samples, with = FALSE])
            rownames(mm) <- counts_sub$gene_name
            
            # Filter genes by row sum
            minsum <- req(input$min_row_sum)
            rs <- rowSums(mm, na.rm = TRUE)
            keep <- which(rs >= minsum)
            if (length(keep) == 0) stop("No genes left after applying min_row_sum filter.")
            mm <- mm[keep, , drop = FALSE]
            
            # Z-score rows if requested
            if (isTRUE(input$scale_rows)) {
                mm <- t(scale(t(mm), center = TRUE, scale = TRUE))
                mm[is.na(mm)] <- 0
            }
            
            # Final meta alignment (in same order as matrix columns)
            ann <- ann[match(colnames(mm), ann$Sample)]
            
            list(mat = mm, ann = ann)
        })
        
        
        build_heatmap <- reactive({
            prep <- heatmap_prep(); req(prep)
            m <- prep$mat; ann <- prep$ann
            
            my_col <- colorRamp2(breaks = c(-4, -2, 0, 2, 4),
                                 colors = c('#00429d', '#73a2c6', 'grey96', '#f4777f', '#93003a'))
            
            # Build the top annotation only if Group exists
            ha <- NULL
            if ("Group" %in% colnames(ann)) {
                uniq <- unique(as.character(ann$Group))
                base_palette <- c("control" = "#3a5cbc", "ibrutinib" = "#a62a17")
                pal <- setNames(
                    ifelse(uniq %in% names(base_palette), base_palette[uniq],
                           grDevices::colorRampPalette(c("#00429d", "#73a2c6", "#f4777f", "#93003a"))(length(uniq))),
                    uniq
                )
                ha <- HeatmapAnnotation(Group = ann$Group, col = list(Group = pal),
                                        annotation_legend_param = list(title = "Group"))
            }
            

            if (input$mode == "supervised") {
                column_split_val <- ann$Group
                cluster_columns_flag <- FALSE
                top_ann <- ha
            } else {
                # Unsupervised: split columns by numeric value if input$col_split >1, else cluster all together
                column_split_val <- as.integer(input$col_split)
                cluster_columns_flag <- TRUE
                # KEEP the annotation so Group shows at the top
                top_ann <- ha
            }
            
            row_split_val <- as.integer(req(input$row_split) %||% 2)
            
            ht <- Heatmap(
                m,
                name = if (isTRUE(input$scale_rows)) "Z-score" else "Value",
                col = my_col,
                use_raster = TRUE,
                clustering_distance_rows = "pearson",
                clustering_distance_columns = "pearson",
                clustering_method_rows = "ward.D2",
                clustering_method_columns = "ward.D2",
                border = TRUE,
                column_split = column_split_val,
                cluster_columns = cluster_columns_flag,
                row_split = if (row_split_val > 1) row_split_val else 1,
                top_annotation = top_ann,
                show_row_names = FALSE
            )
            
            ht
        })
        
        observe({
            req(input$filter_column)
            new_label <- if (input$filter_column == "pvalue") "p-value threshold" else "padj threshold"
            updateNumericInput(session, "pval_thresh", label = new_label)
        })
        
        output$heatmap_plot <- renderPlot({
            ht <- build_heatmap(); req(ht)
            tryCatch({
                grid::grid.newpage()
                ComplexHeatmap::draw(ht, merge_legend = TRUE)
            }, error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("Heatmap error:", e$message), cex = 1)
            })
        }, res = 96)
        

        output$download_heatmap <- downloadHandler(
            filename = function() paste0("heatmap_", Sys.Date(), ".png"),
            content = function(file) {
                ht <- build_heatmap(); req(ht)
                
                # Draw and capture as ggplot
                gr <- grid::grid.grabExpr({
                    ComplexHeatmap::draw(ht, merge_legends = TRUE)
                }) |> ggplotify::as.ggplot()
                
                # Save with fixed size and dpi (wide, high-quality)
                ggplot2::ggsave(
                    filename = file,
                    plot = gr,
                    width = 10, height = 10,
                    units = "in",
                    dpi = 600
                )
            }
        )
        
        
        
    })
}
