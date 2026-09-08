

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

        # Input ID for a covariate's categorical/continuous override toggle.
        # make.names() keeps the ID valid even for column names with spaces/punctuation.
        covtype_input_id <- function(col) paste0("covtype_", make.names(col))

        # Single source of truth for a covariate's fitted type, used identically by the
        # design preview, the analysis run, and the provenance record so they can never
        # disagree. Non-numeric covariates are always categorical (not user-configurable).
        # Numeric covariates default to is_categorical_covariate()'s suggestion but defer
        # to the user's toggle in covariate_type_ui whenever one has been set.
        resolve_covariate_type <- function(col, x) {
            inferred <- if (is.numeric(x)) {
                if (is_categorical_covariate(x)) "categorical" else "continuous"
            } else {
                "categorical"
            }

            override <- if (is.numeric(x)) input[[covtype_input_id(col)]] else NULL
            used <- if (!is.null(override)) override else inferred

            list(
                covariate = col,
                inferred_type = inferred,
                used_type = used,
                was_override = !is.null(override) && !identical(override, inferred)
            )
        }

        # Persists each covariate's categorical/continuous choice across re-renders of
        # covariate_type_ui, keyed by covariate name. covariate_type_ui must not read the
        # toggle inputs directly to decide its own `selected=` value -- that would make the
        # renderUI reactively depend on the very widget it creates, so every click would
        # regenerate the widget back to the inferred default before the click could register.
        # Reading from this store with isolate() breaks that loop while still letting the
        # choice survive when the block regenerates for an unrelated reason (e.g. a second
        # covariate gets added).
        covariate_type_overrides <- reactiveValues()

        observe({
            covs <- input$covariates
            if (is.null(covs)) return()
            for (col in covs) {
                val <- input[[covtype_input_id(col)]]
                if (!is.null(val)) covariate_type_overrides[[col]] <- val
            }
        })

        order_metadata <- function(meta, sample_col, main_factor = NULL) {
            if (!is.null(main_factor) && main_factor %in% colnames(meta)) {
                ord <- order(as.character(meta[[sample_col]]), as.character(meta[[main_factor]]))
            } else {
                ord <- order(as.character(meta[[sample_col]]))
            }
            meta[ord]
        }

        subset_inputs <- reactive({
            meta_df <- meta_data()
            if (is.null(meta_df) ||
                is.null(input$main_factor) ||
                is.null(input$ref_level) ||
                is.null(input$target_level) ||
                identical(input$ref_level, input$target_level)) {
                return(NULL)
            }

            meta <- copy(meta_df)
            samp_col <- get_sample_col(meta)
            keep_levels <- unique(c(input$ref_level, input$target_level))

            subset_meta <- meta[meta[[input$main_factor]] %in% keep_levels, ]
            subset_meta <- subset_meta[!is.na(subset_meta[[samp_col]]), ]

            list(
                sample_col = samp_col,
                keep_levels = keep_levels,
                meta = subset_meta
            )
        })
        
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
                uiOutput(ns("covariate_type_ui")),
                uiOutput(ns("contrast_ui"))
            )

        })

        # The main factor can't also be selected as its own covariate (that produced a
        # duplicated "~ Group + Group" design formula). Whenever the main factor changes,
        # drop it from the covariates choices -- and from the current selection if it was
        # already picked -- without re-rendering the rest of the sidebar.
        observeEvent(input$main_factor, {
            req(meta_data())
            df <- meta_data()
            cols <- colnames(df)
            cols <- cols[!tolower(cols) %in% c("sample", "sampleid", "id", "name")]
            cols <- setdiff(cols, input$main_factor)

            updateSelectizeInput(
                session, "covariates",
                choices = cols,
                selected = setdiff(input$covariates, input$main_factor)
            )
        }, ignoreNULL = TRUE)

        output$covariate_type_ui <- renderUI({

            req(meta_data())
            df <- meta_data()
            ns <- session$ns
            covs <- input$covariates

            if (is.null(covs) || length(covs) == 0) return(NULL)

            rows <- lapply(covs, function(col) {
                if (!col %in% colnames(df)) return(NULL)
                x <- df[[col]]

                if (!is.numeric(x)) {
                    return(tags$div(
                        style = "margin-bottom: 6px; font-size: 0.9em;",
                        strong(col), ": Categorical (fixed, non-numeric)"
                    ))
                }

                inferred <- if (is_categorical_covariate(x)) "categorical" else "continuous"
                n_distinct <- length(unique(x[!is.na(x)]))
                suggested_label <- if (inferred == "categorical") "Categorical" else "Continuous"

                # Restore a previously-made choice for this covariate if there is one;
                # isolate() so this block doesn't get invalidated by the very click it's
                # reacting to (see covariate_type_overrides above).
                current <- isolate(covariate_type_overrides[[col]])
                selected <- if (!is.null(current)) current else inferred

                tags$div(
                    style = "margin-bottom: 10px;",
                    radioButtons(
                        ns(covtype_input_id(col)),
                        label = paste0(col, " (numeric, ", n_distinct, " distinct value",
                                       if (n_distinct != 1) "s" else "", ")"),
                        choices = c("Categorical" = "categorical", "Continuous" = "continuous"),
                        selected = selected,
                        inline = TRUE
                    ),
                    tags$small(style = "color:#888;", paste0("suggested: ", suggested_label))
                )
            })

            tagList(
                hr(),
                h5("Covariate Types"),
                rows
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
                selectInput(ns("target_level"), "Target Level (Treatment):", choices = lvls, selected = lvls[min(2, length(lvls))])
            )
            
        })
        
        
        output$design_preview <- renderText({

            req(input$main_factor)
            factors <- c(input$covariates, input$main_factor)
            preview <- paste("Design Formula:\n~", paste(factors, collapse = " + "))

            df <- meta_data()
            continuous_cols <- character(0)
            overrides <- character(0)

            if (!is.null(df) && length(input$covariates) > 0) {
                for (col in input$covariates) {
                    if (!col %in% colnames(df)) next
                    info <- resolve_covariate_type(col, df[[col]])

                    if (info$used_type == "continuous") continuous_cols <- c(continuous_cols, col)

                    if (info$was_override) {
                        used_label <- if (info$used_type == "categorical") "Categorical" else "Continuous"
                        inferred_label <- if (info$inferred_type == "categorical") "Categorical" else "Continuous"
                        overrides <- c(overrides, paste0(col, " -> ", used_label, " (suggested ", inferred_label, ")"))
                    }
                }
            }

            if (length(continuous_cols) > 0) {
                preview <- paste0(
                    preview,
                    "\n\nFit as continuous (not converted to a factor): ",
                    paste(continuous_cols, collapse = ", ")
                )
            }

            if (length(overrides) > 0) {
                preview <- paste0(preview, "\nOverrides: ", paste(overrides, collapse = "; "))
            }

            preview

        })
        
        
        analysis_out <- eventReactive(input$run_btn, {
            
            req(counts_data(), meta_data(), input$main_factor, input$ref_level, input$target_level)
            
            tryCatch({
                
                cts <- copy(counts_data())
                meta <- copy(meta_data())
                
                gene_col <- colnames(cts)[1] 
                samp_col <- get_sample_col(meta)
                meta <- meta[!is.na(meta[[samp_col]]), ]
                meta <- order_metadata(meta, samp_col, input$main_factor)

                sample_ids <- as.character(meta[[samp_col]])
                overlap <- check_sample_overlap(sample_ids, colnames(cts)[-1],
                                                 meta_label = "metadata", other_label = "counts file")
                if (!is.null(overlap$message)) showNotification(overlap$message, type = "warning", duration = NULL)
                sample_ids <- overlap$common

                meta <- meta[match(sample_ids, as.character(meta[[samp_col]]))]
                keep_cols <- c(gene_col, sample_ids)
                cts <- cts[, ..keep_cols]
                
                mm <- as.matrix(cts[, -1, with=FALSE])
                rownames(mm) <- cts[[1]]
                
                mm <- round(mm)

                validate(
                    need(nrow(meta) >= 2, "Need at least two samples for DESeq2 analysis."),
                    need(!identical(input$ref_level, input$target_level), "Choose two different groups for the comparison."),
                    need(all(c(input$ref_level, input$target_level) %in% unique(as.character(meta[[input$main_factor]]))),
                         "Selected contrast groups were not found in the matched metadata.")
                )
                
                design_cols <- c(input$covariates, input$main_factor)

                # The main factor always needs discrete ref/target levels to relevel()
                # against, so it is not part of the user-configurable type resolution.
                # Covariates are resolved via resolve_covariate_type() -- the same function
                # the design preview uses -- so the fitted model always matches what was
                # shown to the user, including any categorical/continuous override.
                covariate_types <- list()
                for(col in design_cols) {
                    if (col == input$main_factor) {
                        meta[[col]] <- as.factor(meta[[col]])
                    } else {
                        info <- resolve_covariate_type(col, meta[[col]])
                        covariate_types[[col]] <- info
                        if (info$used_type == "categorical") {
                            meta[[col]] <- as.factor(meta[[col]])
                        }
                    }
                }

                design_formula <- as.formula(paste0("~ ", paste(design_cols, collapse = " + ")))
                
                dds <- DESeqDataSetFromMatrix(
                    countData = mm, 
                    colData = meta, 
                    design = design_formula
                )
                
                dds[[input$main_factor]] <- relevel(dds[[input$main_factor]], ref = input$ref_level)
                
                dds <- estimateSizeFactors(dds)
                norm_counts <- counts(dds, normalized = TRUE)
                
                dds <- DESeq(dds)
                
                res <- lfcShrink(
                    dds, 
                    contrast = c(input$main_factor, input$target_level, input$ref_level), 
                    type = "ashr"
                )
                
                res_df <- as.data.frame(res)
                res_df$gene_name <- rownames(res_df)
                res_dt <- as.data.table(res_df)
                setcolorder(res_dt, c("gene_name", setdiff(colnames(res_dt), "gene_name")))
                
                norm_dt <- as.data.frame(norm_counts)
                norm_dt$gene_name <- rownames(norm_dt)
                norm_dt <- as.data.table(norm_dt)
                setcolorder(norm_dt, c("gene_name", setdiff(colnames(norm_dt), "gene_name")))
                
                list(res = res_dt, counts = norm_dt, dds = dds, meta = meta, covariate_types = covariate_types)
                
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
        active_subset <- reactive({
            subset_info <- subset_inputs()
            if (is.null(subset_info) || is.null(counts_data())) return(NULL)

            meta <- copy(subset_info$meta)
            cts <- copy(counts_data())

            gene_col <- colnames(cts)[1]
            keep_sample_ids <- unique(as.character(meta[[subset_info$sample_col]]))
            ordered_sample_ids <- keep_sample_ids[keep_sample_ids %in% colnames(cts)]
            keep_cols <- c(gene_col, ordered_sample_ids)

            if (length(keep_cols) <= 1) return(NULL)

            cts <- cts[, ..keep_cols]

            list(
                counts = cts,
                meta = meta,
                main_factor = input$main_factor,
                ref_level = input$ref_level,
                target_level = input$target_level
            )
        })

        analysis_results <- reactive({
            if (input$run_btn == 0) return(NULL)

            out <- analysis_out()
            if (is.null(out)) return(NULL)

            list(res = out$res, counts = out$counts, meta = out$meta, covariate_types = out$covariate_types)
        })

        return(list(
            results = analysis_results,
            active_subset = active_subset
        ))
        
    })
}

