

deUI <- function(id) {
    ns <- NS(id)
    tagList(
        sidebarLayout(
            sidebarPanel(
                h4("1. Experimental Design"),
                uiOutput(ns("design_controls")),
                hr(),
                verbatimTextOutput(ns("design_preview")),
                actionButton(ns("run_btn"), "Run DESeq2 Analysis", class = "btn-primary", width = "100%"),
                hr(),
                textOutput(ns("error_msg"))
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Results Table",
                             withSpinner(DTOutput(ns("res_table")), type = 4),
                             downloadButton(ns("dl_res"), "Download Results (Excel)")
                    ),
                    tabPanel("Normalized Counts", 
                             withSpinner(DTOutput(ns("norm_table"))),
                             downloadButton(ns("dl_norm"), "Download Normalized Counts")
                    )
                )
            )
        )
    )
}

deserver <- function(id, counts_data, meta_data) {
    moduleServer(id, function(input, output, session) {
        
        output$design_controls <- renderUI({
            req(meta_data()) 
            df <- meta_data()
            ns <- session$ns
            cols <- colnames(df)
            cols <- cols[!tolower(cols) %in% c("sample", "sampleid", "id", "name")]
            tagList(
                selectInput(ns("main_factor"), "Main Variable (Comparison):", choices = cols, selected = "Group"),
                selectizeInput(ns("covariates"), "Covariates (e.g., patientID, Batch):", 
                               choices = cols, multiple = TRUE, options = list(placeholder = "Optional")),
                uiOutput(ns("contrast_ui"))
            )
        })
        
        output$contrast_ui <- renderUI({
            req(input$main_factor, meta_data())
            ns <- session$ns
            df <- meta_data()
            lvls <- unique(as.character(df[[input$main_factor]]))
            tagList(
                hr(),
                h5("Contrast Definition"),
                selectInput(ns("ref_level"), "Reference Level (Control):", choices = lvls, selected = lvls[1]),
                selectInput(ns("target_level"), "Target Level (Treatment):", choices = lvls, selected = lvls[2])
            )
        })
        
        output$design_preview <- renderText({
            req(input$main_factor)
            factors <- c(input$covariates, input$main_factor)
            paste("Design Formula:\n~", paste(factors, collapse = " + "))
        })
        
        analysis_out <- eventReactive(input$run_btn, {
            req(counts_data(), meta_data(), input$main_factor, input$ref_level, input$target_level)
            
            tryCatch({
                cts <- copy(counts_data())
                meta <- copy(meta_data())
                
                gene_col <- colnames(cts)[1] 
                samp_col <- grep("sample|id", colnames(meta), ignore.case = TRUE, value = TRUE)[1]
                if(is.na(samp_col)) samp_col <- colnames(meta)[1]
                
                keep_cols <- c(gene_col, intersect(colnames(cts), meta[[samp_col]]))
                cts <- cts[, ..keep_cols]
                
                mm <- as.matrix(cts[, -1, with=FALSE])
                rownames(mm) <- cts[[1]]
                
                # Rounding to integer (Critical Fix)
                mm <- round(mm)
                
                meta <- meta[meta[[samp_col]] %in% colnames(mm), ]
                meta <- meta[match(colnames(mm), meta[[samp_col]]), ]
                
                design_cols <- c(input$covariates, input$main_factor)
                for(col in design_cols) {
                    meta[[col]] <- as.factor(meta[[col]])
                }
                
                design_formula <- as.formula(paste0("~ ", paste(design_cols, collapse = " + ")))
                
                dds <- DESeqDataSetFromMatrix(countData = mm, colData = meta, design = design_formula)
                dds[[input$main_factor]] <- relevel(dds[[input$main_factor]], ref = input$ref_level)
                
                dds <- estimateSizeFactors(dds)
                norm_counts <- counts(dds, normalized = TRUE)
                
                dds <- DESeq(dds)
                
                res <- lfcShrink(dds, contrast = c(input$main_factor, input$target_level, input$ref_level), type = "ashr")
                
                res_df <- as.data.frame(res)
                res_df$gene_name <- rownames(res_df)
                res_dt <- as.data.table(res_df)
                
                norm_dt <- as.data.frame(norm_counts)
                norm_dt$gene_name <- rownames(norm_dt)
                norm_dt <- as.data.table(norm_dt)
                
                # We return meta as well, so downstream tabs use the aligned metadata
                list(res = res_dt, counts = norm_dt, dds = dds, meta = meta)
                
            }, error = function(e) {
                output$error_msg <- renderText(paste("Analysis Failed:", e$message))
                return(NULL)
            })
        })
        
        output$res_table <- renderDT({
            req(analysis_out())
            datatable(analysis_out()$res[order(padj)], options = list(pageLength = 10, scrollX = TRUE))
        })
        
        output$norm_table <- renderDT({
            req(analysis_out())
            datatable(analysis_out()$counts, options = list(pageLength = 10, scrollX = TRUE))
        })
        
        output$dl_res <- downloadHandler(
            filename = function() { paste0("DESeq2_", input$target_level, "_vs_", input$ref_level, ".csv") },
            content = function(file) { write.csv(analysis_out()$res, file, row.names = FALSE) }
        )
        output$dl_norm <- downloadHandler(
            filename = function() { "gene_counts_deseq2_normalized.csv" },
            content = function(file) { write.csv(analysis_out()$counts, file, row.names = FALSE) }
        )
        
        # --- CRITICAL FIX IS HERE (Step 4) ---
        return(reactive({
            # If the user has NOT clicked the button, return NULL immediately.
            # This tells main.R: "I have no data, go check the uploads instead."
            if (input$run_btn == 0) return(NULL)
            
            # If button clicked, try to get results. 
            # If analysis failed (returned NULL), we pass that NULL along.
            out <- analysis_out()
            if (is.null(out)) return(NULL)
            
            list(res = out$res, counts = out$counts, meta = out$meta)
        }))
    })
}