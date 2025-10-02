# R/heatmap.R
library(shiny)
library(data.table)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(colorRamp2)
library(ggplot2)
library(ggplotify)
# (circlize already used by ComplexHeatmap)

# UI -------------------------------------------------------
heatmapUI <- function(id) {
    ns <- NS(id)
    tagList(
        sidebarLayout(
            sidebarPanel(
                helpText("Heatmap based on DE filtering from DESeq2 (pvalue & log2FoldChange)."),
                numericInput(ns("pval_thresh"), "p-value threshold", value = 0.05, min = 0, step = 0.01),
                numericInput(ns("lfc_thresh"), "abs(log2FoldChange) threshold", value = 1, min = 0, step = 0.1),
                numericInput(ns("min_row_sum"), "Min gene total across samples (rowSums)", value = 10, min = 0, step = 1),
                checkboxInput(ns("scale_rows"), "Scale rows (z-score)", value = TRUE),
                numericInput(ns("col_split"), "column_split (default 2)", value = 2, min = 1, step = 1),
                numericInput(ns("row_split"), "row_split (default 2)", value = 2, min = 1, step = 1),
                hr(),
                downloadButton(ns("download_heatmap"), "Download Heatmap (PNG)")
            ),
            mainPanel(
                plotOutput(ns("heatmap_plot"), height = "800px")
            )
        )
    )
}

# SERVER ---------------------------------------------------
heatmapServer <- function(id, meta_data, counts_data, deseq_data) {
    moduleServer(id, function(input, output, session) {
        
        # Reactive copies
        hm_data <- reactive({
            req(meta_data(), counts_data(), deseq_data())
            meta <- copy(meta_data())
            counts <- copy(counts_data())
            deseq <- copy(deseq_data())
            
            # Tolerant renaming for common variants:
            if ("Geneid" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) {
                data.table::setnames(counts, "Geneid", "gene_name")
            }
            if ("gene_id" %in% colnames(counts) && !"gene_name" %in% colnames(counts)) {
                data.table::setnames(counts, "gene_id", "gene_name")
            }
            # metadata sample column name: prefer 'Sample' but accept sampleID
            if (!("Sample" %in% colnames(meta)) && ("sampleID" %in% colnames(meta))) {
                data.table::setnames(meta, "sampleID", "Sample")
            }
            
            list(meta = meta, counts = counts, deseq = deseq)
        })
        
        # Prepare matrix following original script logic
        heatmap_prep <- reactive({
            dd <- hm_data()
            req(dd)
            meta <- dd$meta
            counts <- dd$counts
            deseq <- dd$deseq
            
            # Validate presence of necessary columns
            req("Sample" %in% colnames(meta), "metadata must contain a 'Sample' column")
            req("gene_name" %in% colnames(counts), "counts must contain a 'gene_name' column")
            req(("pvalue" %in% colnames(deseq) && "log2FoldChange" %in% colnames(deseq)),
                "deseq table must contain 'pvalue' and 'log2FoldChange'")
            
            # Filter DE genes by thresholds (defaults match your script)
            pth <- req(input$pval_thresh); lfc <- req(input$lfc_thresh)
            d_sig <- deseq[which(pvalue <= pth & abs(log2FoldChange) >= lfc)]
            gene_index <- sort(unique(d_sig$GeneID %||% d_sig$Geneid %||% d_sig$gene_id))
            if (length(gene_index) == 0) stop("No significant genes found with the chosen thresholds.")
            
            # Subset counts to the selected genes
            counts_sub <- counts[gene_name %in% gene_index]
            if (nrow(counts_sub) == 0) stop("No matching genes found in counts for DE gene list.")
            
            # Select metadata rows corresponding to columns in counts_sub
            sample_cols <- setdiff(colnames(counts_sub), "gene_name")
            ann <- meta[Sample %in% sample_cols]
            
            # Set Group factor levels as in script (if Group exists), else leave as-is
            if ("Group" %in% colnames(ann)) {
                # if control & ibrutinib expected, keep that order; otherwise preserve factor levels if present
                if (all(c("control", "ibrutinib") %in% unique(as.character(ann$Group)))) {
                    ann$Group <- factor(ann$Group, levels = c("control", "ibrutinib"))
                } else {
                    ann$Group <- as.factor(ann$Group)
                }
                ann <- ann[order(ann$Group)]
            } else {
                ann <- ann[order(ann$Sample)]
            }
            
            # Reorder columns of counts accordingly
            sel_samples <- intersect(ann$Sample, sample_cols)
            if (length(sel_samples) == 0) stop("No sample names overlap between metadata$Sample and counts columns.")
            mm <- as.matrix(counts_sub[, sel_samples, with = FALSE])
            rownames(mm) <- counts_sub$gene_name
            
            # Filter genes with low total counts
            minsum <- req(input$min_row_sum)
            rs <- rowSums(mm, na.rm = TRUE)
            keep <- which(rs >= minsum)
            if (length(keep) == 0) stop("No genes left after applying min_row_sum filter.")
            mm <- mm[keep, , drop = FALSE]
            
            # Scale rows (z-score) if requested
            if (isTRUE(input$scale_rows)) {
                mm <- t(scale(t(mm), center = TRUE, scale = TRUE))
                mm[is.na(mm)] <- 0
            }
            
            # Ensure meta matches final columns order
            ann <- ann[match(colnames(mm), ann$Sample)]
            list(mat = mm, ann = ann)
        })
        
        # Build ComplexHeatmap object
        build_heatmap <- reactive({
            prep <- heatmap_prep()
            req(prep)
            m <- prep$mat
            ann <- prep$ann
            
            # color mapping from your script
            my_col <- colorRamp2(breaks = c(-4, -2, 0, 2, 4),
                                 colors = c('#00429d', '#73a2c6', 'grey96', '#f4777f', '#93003a'))
            
            # Column annotation based on Group (if present)
            if ("Group" %in% colnames(ann)) {
                # preserve control/ibrutinib colors if present; else generate palette
                uniq <- unique(as.character(ann$Group))
                base_palette <- c("control" = "#3a5cbc", "ibrutinib" = "#a62a17")
                pal <- setNames(
                    ifelse(uniq %in% names(base_palette), base_palette[uniq],
                           grDevices::colorRampPalette(c("#00429d", "#73a2c6", "#f4777f", "#93003a"))(length(uniq))),
                    uniq
                )
                ha <- HeatmapAnnotation(Group = ann$Group, col = list(Group = pal),
                                        annotation_legend_param = list(title = "Group"))
            } else {
                ha <- HeatmapAnnotation(dummy = anno_empty())
            }
            
            col_split_val <- as.integer(req(input$col_split) %||% 2)
            row_split_val <- as.integer(req(input$row_split) %||% 2)
            
            ht <- Heatmap(
                m,
                name = if (isTRUE(input$scale_rows)) "Z-score" else "Value",
                col = my_col,
                use_raster = TRUE,
                clustering_distance_rows = "euclidean",
                clustering_distance_columns = "euclidean",
                clustering_method_rows = "ward.D2",
                clustering_method_columns = "ward.D2",
                border = TRUE,
                column_split = if (col_split_val > 1) col_split_val else 1,
                row_split = if (row_split_val > 1) row_split_val else 1,
                top_annotation = ha,
                show_row_names = FALSE
            )
            
            ht
        })
        
        # Render heatmap
        output$heatmap_plot <- renderPlot({
            ht <- build_heatmap()
            req(ht)
            tryCatch({
                grid::grid.newpage()
                ComplexHeatmap::draw(ht, merge_legend = TRUE)
            }, error = function(e) {
                plot.new()
                text(0.5, 0.5, paste("Heatmap error:", e$message), cex = 1)
            })
        }, res = 96)
        
        # Download handler (PNG)
        output$download_heatmap <- downloadHandler(
            filename = function() {
                paste0("heatmap_", Sys.Date(), ".png")
            },
            content = function(file) {
                ht <- build_heatmap()
                req(ht)
                prep <- heatmap_prep()
                nr <- nrow(prep$mat); nc <- ncol(prep$mat)
                # simple heuristics for PNG size
                png_w <- min(4000, max(1200, nc * 40))
                png_h <- min(4000, max(1200, nr * 6))
                png(filename = file, width = png_w, height = png_h, res = 300)
                tryCatch({
                    ComplexHeatmap::draw(ht, merge_legend = TRUE)
                }, finally = {
                    dev.off()
                })
            }
        )
        
    })
}
