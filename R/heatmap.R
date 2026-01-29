
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
                
                h4("Visual Parameters"),
                # New parameters to control visuals like in your script
                numericInput(ns("font_size"), "Font Size", value = 12, min = 6),
                checkboxInput(ns("show_row_names"), "Show Row Names", value = FALSE),
                
                conditionalPanel(condition = sprintf("input['%s'] == 'unsupervised'", ns("mode")),
                                 numericInput(ns("col_split"), "Column Split (k)", value = 2, min = 1, step = 1)
                ),
                numericInput(ns("row_split"), "Row Split (k)", value = 2, min = 1, step = 1),
                
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
            
            # Tolerant renaming
            if ("Geneid" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) data.table::setnames(counts, "Geneid", "gene_name")
            if ("gene_id" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) data.table::setnames(counts, "gene_id", "gene_name")
            if (!("Sample" %in% colnames(meta)) && ("sampleID" %in% colnames(meta))) data.table::setnames(meta, "sampleID", "Sample")
            
            list(meta = meta, counts = counts, deseq = deseq)
        })
        
        # --- 2. Data Preparation (The "Global" Logic) ---
        heatmap_prep <- reactive({
            dd <- hm_data(); req(dd)
            meta <- dd$meta; counts <- dd$counts; deseq <- dd$deseq
            
            # A. Filter Significant Genes
            # (We do this BEFORE filtering samples to match script)
            filter_col <- req(input$filter_column)
            
            # FIX: Use strict '>' to match script logic exactly
            d_sig <- deseq[get(filter_col) <= input$pval_thresh & abs(log2FoldChange) > input$lfc_thresh]
            
            gene_index <- sort(unique(d_sig$GeneID %||% d_sig$Geneid %||% d_sig$gene_name))
            if (length(gene_index) == 0) stop("No significant genes found.")
            
            # B. Filter Counts (Keep ALL samples for now)
            counts_filtered <- counts[gene_name %in% gene_index]
            if (nrow(counts_filtered) == 0) stop("No matching genes in counts.")
            
            # C. Prepare Matrix
            mm <- as.matrix(counts_filtered[, -1]) # remove gene_name
            rownames(mm) <- counts_filtered$gene_name
            
            # D. Filter Low Counts (Global)
            # Applied to all samples before subsetting
            keep_rows <- which(rowSums(mm, na.rm=TRUE) >= input$min_row_sum)
            mm <- mm[keep_rows, , drop = FALSE]
            if (nrow(mm) == 0) stop("All genes filtered by row sum.")
            
            # E. Calculate Z-Score Globally (The Script Logic)
            # We scale first, THEN subset. This makes colors comparable across comparisons.
            if (isTRUE(input$scale_rows)) {
                mm <- t(scale(t(mm), center = TRUE, scale = TRUE))
                mm[is.na(mm)] <- 0
            }
            
            # F. Now Subset Samples for the View
            # 1. Identify valid samples from Metadata
            keep_samples <- intersect(meta$Sample, colnames(mm))
            if (length(keep_samples) == 0) stop("No overlapping samples found.")
            
            # 2. Subset the ALREADY SCALED matrix
            mm_view <- mm[, keep_samples]
            
            # 3. Align Metadata
            ann <- meta[Sample %in% keep_samples]
            ann <- ann[match(colnames(mm_view), ann$Sample)]
            
            # If Supervised, order by Group
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
            
            # FIX: Set Seed for Consistent Unsupervised Clustering
            set.seed(123)
            
            # Colors (Matches script)
            my_col <- colorRamp2(breaks = c(-4, -2, 0, 2, 4),
                                 colors = c('#00429d', '#73a2c6', 'grey96', '#f4777f', '#93003a'))
            
            # Annotation
            # ha <- NULL
            # if ("Group" %in% colnames(ann)) {
            #     uniq <- unique(as.character(ann$Group))
            #     base_palette <- c("control" = "#3a5cbc", "ibrutinib" = "#a62a17",
            #                       "Subsets_1_99" = "#a62a17", "Other_U_CLL" = "#3a5cbc")
            #     
            #     pal <- setNames(
            #         ifelse(uniq %in% names(base_palette), base_palette[uniq],
            #                grDevices::colorRampPalette(c("#00429d", "#73a2c6", "#f4777f", "#93003a"))(length(uniq))),
            #         uniq
            #     )
            #     ha <- HeatmapAnnotation(Group = ann$Group, col = list(Group = pal),
            #                             simple_anno_size = unit(0.3, "cm"),
            #                             annotation_legend_param = list(title = "Group"))
            # }
            
            # --- Dynamic Annotation Colors ---
            ha <- NULL
            if ("Group" %in% colnames(ann)) {
                
                # 1. Get all unique groups in the correct order (Factor levels)
                # If it's already a factor, levels() gives the correct order (Ref first).
                # If it's character, unique() gives order of appearance.
                if (is.factor(ann$Group)) {
                    groups <- levels(ann$Group)
                } else {
                    groups <- unique(as.character(ann$Group))
                }
                
                # 2. Define your two main custom colors
                col_ref <- "#3a5cbc" # Blue (Control)
                col_trt <- "#a62a17" # Red (Treatment)
                
                # 3. Assign colors dynamically
                if (length(groups) == 2) {
                    # Exactly 2 groups: First gets Blue, Second gets Red
                    final_colors <- c(col_ref, col_trt)
                } else {
                    # More (or less) than 2 groups: Create a gradient/palette 
                    # that spans from Blue to Red (or pick distinct colors)
                    final_colors <- colorRampPalette(c(col_ref, "#f4777f", col_trt))(length(groups))
                }
                
                # 4. Map the colors to the names
                # This creates a named vector: c("GroupA" = Blue, "GroupB" = Red)
                pal <- setNames(final_colors, groups)
                
                # 5. Create Annotation
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
            
            # Plot
            Heatmap(
                m,
                name = if (isTRUE(input$scale_rows)) "Z-score" else "Value",
                col = my_col,
                use_raster = TRUE,
                
                clustering_distance_rows = "pearson",
                clustering_distance_columns = "pearson",
                clustering_method_rows = "ward.D2",
                clustering_method_columns = "ward.D2",
                
                border = TRUE,
                rect_gp = gpar(col = "white", lwd = .25),
                
                column_split = col_split_val,
                cluster_columns = cluster_cols,
                row_split = row_split_val,
                
                top_annotation = ha,
                
                # Visuals
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
                # Ensure seed is set for download too
                set.seed(123)
                ht <- build_heatmap()
                gr <- grid::grid.grabExpr(ComplexHeatmap::draw(ht, merge_legends = TRUE)) |> ggplotify::as.ggplot()
                ggplot2::ggsave(file, plot = gr, width = 10, height = 10, dpi = 600)
            }
        )
        
        observe({
            req(input$filter_column)
            updateNumericInput(session, "pval_thresh", label = if (input$filter_column == "pvalue") "p-value threshold" else "padj threshold")
        })
    })
}
# 
# heatmapUI <- function(id) {
#     ns <- NS(id)
#     tagList(
#         sidebarLayout(
#             sidebarPanel(
#                 radioButtons(ns("mode"), "Mode:",
#                              choices = c("Supervised (by Group)" = "supervised",
#                                          "Unsupervised (cluster columns)" = "unsupervised"),
#                              selected = "supervised"),
#                 radioButtons(
#                     ns("filter_column"),
#                     "Filter DE genes by:",
#                     choices = c("pvalue" = "pvalue", "padj" = "padj"),
#                     selected = "pvalue"
#                 ),
#                 numericInput(ns("pval_thresh"), "p-value threshold", value = 0.05, min = 0, step = 0.01),
#                 numericInput(ns("lfc_thresh"), "abs(log2FoldChange) threshold", value = 1, min = 0, step = 0.1),
#                 # numericInput(ns("min_row_sum"), "Min gene total across samples (rowSums)", value = 10, min = 0, step = 1),
#                 checkboxInput(ns("scale_rows"), "Scale rows (z-score)", value = TRUE),
#                 conditionalPanel(condition = sprintf("input['%s'] == 'unsupervised'", ns("mode")),
#                                  numericInput(ns("col_split"), "column_split (numeric)", value = 2, min = 1, step = 1)
#                 ),
#                 numericInput(ns("row_split"), "Row split", value = 2, min = 1, step = 1),
#                 hr(),
#                 downloadButton(ns("download_heatmap"), "Download Heatmap (PNG)")
#             ),
#             mainPanel(
#                 plotOutput(ns("heatmap_plot"), height = "800px")
#             )
#         )
#     )
# }
# 
# heatmapServer <- function(id, meta_data, counts_data, deseq_data) {
#     moduleServer(id, function(input, output, session) {
# 
#         hm_data <- reactive({
#             req(meta_data(), counts_data(), deseq_data())
#             meta <- copy(meta_data())
#             counts <- copy(counts_data())
#             deseq <- copy(deseq_data())
# 
#             # tolerant renaming
#             if ("Geneid" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) {
#                 data.table::setnames(counts, "Geneid", "gene_name")
#             }
#             if ("gene_id" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) {
#                 data.table::setnames(counts, "gene_id", "gene_name")
#             }
#             if (!("Sample" %in% colnames(meta)) && ("sampleID" %in% colnames(meta))) {
#                 data.table::setnames(meta, "sampleID", "Sample")
#             }
# 
#             list(meta = meta, counts = counts, deseq = deseq)
#         })
# 
#         heatmap_prep <- reactive({
#             dd <- hm_data(); req(dd)
#             meta <- dd$meta; counts <- dd$counts; deseq <- dd$deseq
# 
#             # validate
#             req("Sample" %in% colnames(meta), "metadata must contain a 'Sample' column")
#             req("gene_name" %in% colnames(counts), "counts must contain 'gene_name' column")
#             req(("pvalue" %in% colnames(deseq) && "log2FoldChange" %in% colnames(deseq)),
#                 "deseq must contain 'pvalue' and 'log2FoldChange' columns")
# 
#             # Determine which column to filter on (pvalue or padj)
#             filter_col <- req(input$filter_column)  # user choice: "pvalue" or "padj"
#             if (!(filter_col %in% colnames(deseq))) stop(paste0("Column '", filter_col, "' not found in DESeq2 results."))
# 
#             # Filter DE genes
#             pth <- req(input$pval_thresh)
#             lfc <- req(input$lfc_thresh)
#             d_sig <- deseq[get(filter_col) <= pth & abs(log2FoldChange) >= lfc]
# 
#             gene_index <- sort(unique(d_sig$GeneID %||% d_sig$Geneid %||% d_sig$gene_id))
#             if (length(gene_index) == 0) stop("No significant genes found with the chosen thresholds.")
# 
#             counts_sub <- counts[gene_name %in% gene_index]
#             if (nrow(counts_sub) == 0) stop("No matching genes found in counts for DE gene list.")
# 
#             # Select metadata rows corresponding to columns in counts_sub
#             sample_cols <- setdiff(colnames(counts_sub), "gene_name")
#             ann <- meta[Sample %in% sample_cols]
# 
#             # For supervised, ensure Group exists and set factor levels like your script
#             if (input$mode == "supervised") {
#                 if (!("Group" %in% colnames(ann))) stop("Supervised mode requires 'Group' column in metadata.")
#                 # enforce factor levels if both known
#                 if (all(c("control", "ibrutinib") %in% unique(as.character(ann$Group)))) {
#                     ann$Group <- factor(ann$Group, levels = c("control", "ibrutinib"))
#                 } else {
#                     ann$Group <- as.factor(ann$Group)
#                 }
#                 # order samples by Group
#                 ann <- ann[order(ann$Group)]
#             } else {
#                 # Unsupervised: keep metadata ordering as-is (or by Sample)
#                 ann <- ann[order(ann$Sample)]
#             }
# 
#             # Reorder counts columns to match ann$Sample
#             sel_samples <- intersect(ann$Sample, sample_cols)
#             if (length(sel_samples) == 0) stop("No sample names overlap between metadata$Sample and counts columns.")
#             mm <- as.matrix(counts_sub[, sel_samples, with = FALSE])
#             rownames(mm) <- counts_sub$gene_name
# 
#             # Filter genes by row sum
#             # minsum <- req(input$min_row_sum)
#             rs <- rowSums(mm, na.rm = TRUE)
#             # keep <- which(rs >= minsum)
#             # if (length(keep) == 0) stop("No genes left after applying min_row_sum filter.")
#             # mm <- mm[keep, , drop = FALSE]
# 
#             # Z-score rows if requested
#             if (isTRUE(input$scale_rows)) {
#                 mm <- t(scale(t(mm), center = TRUE, scale = TRUE))
#                 mm[is.na(mm)] <- 0
#             }
# 
#             # Final meta alignment (in same order as matrix columns)
#             ann <- ann[match(colnames(mm), ann$Sample)]
# 
#             list(mat = mm, ann = ann)
#         })
# 
# 
#         build_heatmap <- reactive({
#             prep <- heatmap_prep(); req(prep)
#             m <- prep$mat; ann <- prep$ann
# 
#             my_col <- colorRamp2(breaks = c(-4, -2, 0, 2, 4),
#                                  colors = c('#00429d', '#73a2c6', 'grey96', '#f4777f', '#93003a'))
# 
#             # Build the top annotation only if Group exists
#             ha <- NULL
#             if ("Group" %in% colnames(ann)) {
#                 uniq <- unique(as.character(ann$Group))
#                 base_palette <- c("control" = "#3a5cbc", "ibrutinib" = "#a62a17")
#                 pal <- setNames(
#                     ifelse(uniq %in% names(base_palette), base_palette[uniq],
#                            grDevices::colorRampPalette(c("#00429d", "#73a2c6", "#f4777f", "#93003a"))(length(uniq))),
#                     uniq
#                 )
#                 ha <- HeatmapAnnotation(Group = ann$Group, col = list(Group = pal),
#                                         annotation_legend_param = list(title = "Group"))
#             }
# 
# 
#             if (input$mode == "supervised") {
#                 column_split_val <- ann$Group
#                 cluster_columns_flag <- FALSE
#                 top_ann <- ha
#             } else {
#                 # Unsupervised: split columns by numeric value if input$col_split >1, else cluster all together
#                 column_split_val <- as.integer(input$col_split)
#                 cluster_columns_flag <- TRUE
#                 # KEEP the annotation so Group shows at the top
#                 top_ann <- ha
#             }
# 
#             row_split_val <- as.integer(req(input$row_split) %||% 2)
# 
#             ht <- Heatmap(
#                 m,
#                 name = if (isTRUE(input$scale_rows)) "Z-score" else "Value",
#                 col = my_col,
#                 use_raster = TRUE,
#                 clustering_distance_rows = "pearson",
#                 clustering_distance_columns = "pearson",
#                 clustering_method_rows = "ward.D2",
#                 clustering_method_columns = "ward.D2",
#                 border = TRUE,
#                 column_split = column_split_val,
#                 cluster_columns = cluster_columns_flag,
#                 row_split = if (row_split_val > 1) row_split_val else 1,
#                 top_annotation = top_ann,
#                 show_row_names = FALSE
#             )
# 
#             ht
#         })
# 
#         observe({
#             req(input$filter_column)
#             new_label <- if (input$filter_column == "pvalue") "p-value threshold" else "padj threshold"
#             updateNumericInput(session, "pval_thresh", label = new_label)
#         })
# 
#         output$heatmap_plot <- renderPlot({
#             ht <- build_heatmap(); req(ht)
#             tryCatch({
#                 grid::grid.newpage()
#                 ComplexHeatmap::draw(ht, merge_legend = TRUE)
#             }, error = function(e) {
#                 plot.new()
#                 text(0.5, 0.5, paste("Heatmap error:", e$message), cex = 1)
#             })
#         }, res = 96)
# 
# 
#         output$download_heatmap <- downloadHandler(
#             filename = function() paste0("heatmap_", Sys.Date(), ".png"),
#             content = function(file) {
#                 ht <- build_heatmap(); req(ht)
# 
#                 # Draw and capture as ggplot
#                 gr <- grid::grid.grabExpr({
#                     ComplexHeatmap::draw(ht, merge_legends = TRUE)
#                 }) |> ggplotify::as.ggplot()
# 
#                 # Save with fixed size and dpi (wide, high-quality)
#                 ggplot2::ggsave(
#                     filename = file,
#                     plot = gr,
#                     width = 10, height = 10,
#                     units = "in",
#                     dpi = 600
#                 )
#             }
#         )
# 
# 
# 
#     })
# }
# 
# heatmapUI <- function(id) {
#     ns <- NS(id)
#     tagList(
#         sidebarLayout(
#             sidebarPanel(
#                 # --- Section 1: Data Logic ---
#                 h4("Data Parameters"),
#                 radioButtons(ns("mode"), "Mode:",
#                              choices = c("Supervised (by Group)" = "supervised",
#                                          "Unsupervised (cluster columns)" = "unsupervised"),
#                              selected = "supervised"),
#                 radioButtons(ns("filter_column"), "Filter DE genes by:",
#                              choices = c("pvalue" = "pvalue", "padj" = "padj"), selected = "pvalue"),
#                 numericInput(ns("pval_thresh"), "p-value threshold", value = 0.05, min = 0, step = 0.01),
#                 numericInput(ns("lfc_thresh"), "abs(log2FoldChange) threshold", value = 1, min = 0, step = 0.1),
#                 numericInput(ns("min_row_sum"), "Min gene total (rowSums)", value = 10, min = 0, step = 1),
#                 checkboxInput(ns("scale_rows"), "Scale rows (z-score)", value = TRUE),
# 
#                 hr(),
# 
#                 # --- Section 2: Visual Parameters (NEW) ---
#                 h4("Visual Parameters"),
# 
#                 # Clustering options
#                 selectInput(ns("clust_method"), "Clustering Method",
#                             choices = c("ward.D2", "complete", "average", "single"), selected = "ward.D2"),
#                 selectInput(ns("dist_method"), "Distance Method",
#                             choices = c("pearson", "euclidean", "manhattan"), selected = "pearson"),
# 
#                 # Split options
#                 conditionalPanel(condition = sprintf("input['%s'] == 'unsupervised'", ns("mode")),
#                                  numericInput(ns("col_split"), "Column Split (k)", value = 2, min = 1, step = 1)
#                 ),
#                 numericInput(ns("row_split"), "Row Split (k)", value = 2, min = 1, step = 1),
# 
#                 # Font & Label options
#                 checkboxInput(ns("show_row_names"), "Show Row Names", value = FALSE),
#                 numericInput(ns("font_size"), "Font Size (Base)", value = 12, min = 6, max = 24),
# 
#                 hr(),
#                 downloadButton(ns("download_heatmap"), "Download Heatmap (PNG)")
#             ),
#             mainPanel(
#                 plotOutput(ns("heatmap_plot"), height = "800px")
#             )
#         )
#     )
# }
# 
# 
# 
# heatmapServer <- function(id, meta_data, counts_data, deseq_data) {
#     moduleServer(id, function(input, output, session) {
# 
#         # --- 1. Data Cleaning ---
#         hm_data <- reactive({
#             req(meta_data(), counts_data(), deseq_data())
#             meta <- copy(meta_data())
#             counts <- copy(counts_data())
#             deseq <- copy(deseq_data())
# 
#             # Tolerant renaming
#             if ("Geneid" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) data.table::setnames(counts, "Geneid", "gene_name")
#             if ("gene_id" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) data.table::setnames(counts, "gene_id", "gene_name")
#             if (!("Sample" %in% colnames(meta)) && ("sampleID" %in% colnames(meta))) data.table::setnames(meta, "sampleID", "Sample")
# 
#             list(meta = meta, counts = counts, deseq = deseq)
#         })
# 
#         # --- 2. Data Prep (Filtering & Scaling) ---
#         heatmap_prep <- reactive({
#             dd <- hm_data(); req(dd)
#             meta <- dd$meta; counts <- dd$counts; deseq <- dd$deseq
# 
#             # Sync Counts/Meta
#             keep_samples <- intersect(meta$Sample, colnames(counts))
#             if (length(keep_samples) == 0) stop("No overlapping samples found.")
# 
#             counts_sub <- counts[, c("gene_name", keep_samples), with = FALSE]
#             ann <- meta[Sample %in% keep_samples]
# 
#             # Filter Genes
#             filter_col <- req(input$filter_column)
#             if (!(filter_col %in% colnames(deseq))) stop("Filter column missing in DESeq.")
# 
#             d_sig <- deseq[get(filter_col) <= input$pval_thresh & abs(log2FoldChange) >= input$lfc_thresh]
#             gene_index <- sort(unique(d_sig$GeneID %||% d_sig$Geneid %||% d_sig$gene_name))
# 
#             if (length(gene_index) == 0) stop("No significant genes found.")
#             counts_sub <- counts_sub[gene_name %in% gene_index]
#             if (nrow(counts_sub) == 0) stop("No matching genes in counts.")
# 
#             # Matrix Prep
#             mm <- as.matrix(counts_sub[, -1])
#             rownames(mm) <- counts_sub$gene_name
# 
#             # Min Row Sum
#             keep_rows <- which(rowSums(mm, na.rm=TRUE) >= input$min_row_sum)
#             mm <- mm[keep_rows, , drop = FALSE]
#             if (nrow(mm) == 0) stop("All genes filtered by row sum.")
# 
#             # Z-Score (Local)
#             if (isTRUE(input$scale_rows)) {
#                 mm <- t(scale(t(mm), center = TRUE, scale = TRUE))
#                 mm[is.na(mm)] <- 0
#             }
# 
#             # Align
#             ann <- ann[match(colnames(mm), ann$Sample)]
#             if (input$mode == "supervised" && "Group" %in% colnames(ann)) {
#                 ann <- ann[order(ann$Group)]
#                 mm <- mm[, ann$Sample]
#             }
# 
#             list(mat = mm, ann = ann)
#         })
# 
#         # --- 3. Build ComplexHeatmap ---
#         build_heatmap <- reactive({
#             prep <- heatmap_prep(); req(prep)
#             m <- prep$mat
#             ann <- prep$ann
# 
#             # --- YOUR CUSTOM COLORS (PRESERVED) ---
#             my_col <- colorRamp2(breaks = c(-4, -2, 0, 2, 4),
#                                  colors = c('#00429d', '#73a2c6', 'grey96', '#f4777f', '#93003a'))
# 
#             # Top Annotation
#             ha <- NULL
#             if ("Group" %in% colnames(ann)) {
#                 uniq <- unique(as.character(ann$Group))
#                 base_palette <- c("control" = "#3a5cbc", "ibrutinib" = "#a62a17",
#                                   "Subsets_1_99" = "#a62a17", "Other_U_CLL" = "#3a5cbc")
# 
#                 pal <- setNames(
#                     ifelse(uniq %in% names(base_palette), base_palette[uniq],
#                            grDevices::colorRampPalette(c("#00429d", "#73a2c6", "#f4777f", "#93003a"))(length(uniq))),
#                     uniq
#                 )
#                 ha <- HeatmapAnnotation(Group = ann$Group, col = list(Group = pal),
#                                         simple_anno_size = unit(0.3, "cm"),
#                                         annotation_legend_param = list(title = "Group"))
#             }
# 
#             # Logic for Splitting
#             if (input$mode == "supervised") {
#                 col_split_val <- if("Group" %in% colnames(ann)) ann$Group else NULL
#                 cluster_cols <- FALSE
#             } else {
#                 col_split_val <- as.integer(input$col_split)
#                 cluster_cols <- TRUE
#             }
# 
#             row_split_val <- as.integer(input$row_split %||% 1)
#             if(row_split_val < 2) row_split_val <- NULL
# 
#             # --- DRAW HEATMAP ---
#             Heatmap(
#                 m,
#                 name = if (isTRUE(input$scale_rows)) "Z-score" else "Value",
#                 col = my_col,
#                 use_raster = TRUE,
# 
#                 # User Parameters for Clustering
#                 clustering_distance_rows = input$dist_method,
#                 clustering_distance_columns = input$dist_method,
#                 clustering_method_rows = input$clust_method,
#                 clustering_method_columns = input$clust_method,
# 
#                 border = TRUE,
#                 rect_gp = gpar(col = "white", lwd = .25),
# 
#                 column_split = col_split_val,
#                 cluster_columns = cluster_cols,
#                 row_split = row_split_val,
# 
#                 top_annotation = ha,
# 
#                 # User Parameters for Visuals
#                 show_row_names = input$show_row_names,
#                 row_names_gp = gpar(fontsize = input$font_size * 0.8),
#                 column_names_gp = gpar(fontsize = input$font_size * 0.8),
#                 heatmap_legend_param = list(labels_gp = gpar(fontsize = input$font_size))
#             )
#         })
# 
#         # --- 4. Render & Download ---
#         output$heatmap_plot <- renderPlot({
#             ht <- build_heatmap(); req(ht)
#             grid::grid.newpage()
#             ComplexHeatmap::draw(ht, merge_legend = TRUE)
#         }, res = 96)
# 
#         observe({
#             req(input$filter_column)
#             updateNumericInput(session, "pval_thresh", label = if (input$filter_column == "pvalue") "p-value threshold" else "padj threshold")
#         })
# 
#         output$download_heatmap <- downloadHandler(
#             filename = function() paste0("heatmap_", Sys.Date(), ".png"),
#             content = function(file) {
#                 ht <- build_heatmap()
#                 gr <- grid::grid.grabExpr(ComplexHeatmap::draw(ht, merge_legends = TRUE)) |> ggplotify::as.ggplot()
#                 ggplot2::ggsave(file, plot = gr, width = 10, height = 10, dpi = 600)
#             }
#         )
#     })
# }