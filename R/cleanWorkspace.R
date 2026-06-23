# Launch a Shiny app to interactively remove large objects from the workspace.
# Call cleanWorkspace() to open the app. Environment is rebound to shiny's
# namespace at load time (see R/zzz.R) so the unqualified shiny UI/server
# functions resolve without attaching shiny.

cleanWorkspace <- function(env = .GlobalEnv) {

    format_size <- function(bytes) {
        if (is.na(bytes)) return("?")
        if (bytes >= 1e9) sprintf("%.2f GB", bytes / 1e9)
        else if (bytes >= 1e6) sprintf("%.2f MB", bytes / 1e6)
        else if (bytes >= 1e3) sprintf("%.1f KB", bytes / 1e3)
        else sprintf("%d B", as.integer(bytes))
    }

    get_obj_table <- function() {
        nms <- ls(envir = env)
        if (length(nms) == 0) return(data.frame(name = character(), size_bytes = numeric(),
                                                 size_fmt = character(), class = character(),
                                                 stringsAsFactors = FALSE))
        sizes <- sapply(nms, function(nm) {
            tryCatch(as.numeric(object.size(get(nm, envir = env))), error = function(e) NA_real_)
        })
        classes <- sapply(nms, function(nm) {
            tryCatch(class(get(nm, envir = env))[1], error = function(e) "?")
        })
        df <- data.frame(name = nms, size_bytes = sizes, class = classes,
                         stringsAsFactors = FALSE)
        df <- df[order(df$size_bytes, decreasing = TRUE), ]
        df$size_fmt <- sapply(df$size_bytes, format_size)
        rownames(df) <- NULL
        df
    }

    ui <- fluidPage(
        tags$head(tags$style(HTML("
            .obj-row { padding: 4px 0; border-bottom: 1px solid #eee; }
            .obj-row:hover { background-color: #f8f9fa; }
            .header-row { font-weight: bold; border-bottom: 2px solid #dee2e6; padding: 6px 0; }
            .size-col { color: #555; font-family: monospace; }
            .class-col { color: #888; font-style: italic; }
            .btn-danger { margin-top: 4px; }
            body { padding: 20px; }
        "))),
        h3("Workspace Cleaner"),
        fluidRow(
            column(12,
                actionButton("refresh",    "Refresh",          icon = icon("rotate"),  class = "btn-secondary"),
                actionButton("select_all", "Select All",       icon = icon("check-square")),
                actionButton("desel_all",  "Deselect All",     icon = icon("square")),
                actionButton("remove_btn", "Remove Selected",  icon = icon("trash"),   class = "btn-danger"),
                hr()
            )
        ),
        fluidRow(
            column(12, uiOutput("obj_list"))
        )
    )

    server <- function(input, output, session) {

        obj_data <- reactiveVal(get_obj_table())

        observeEvent(input$refresh, {
            obj_data(get_obj_table())
        })

        output$obj_list <- renderUI({
            df <- obj_data()

            if (nrow(df) == 0) return(p("Environment is empty.", style = "color: grey;"))

            header <- fluidRow(class = "header-row",
                column(1, ""),
                column(4, "Name"),
                column(2, "Size"),
                column(2, "Class")
            )

            rows <- lapply(seq_len(nrow(df)), function(i) {
                nm <- df$name[i]
                fluidRow(class = "obj-row",
                    column(1, checkboxInput(paste0("chk_", nm), label = NULL, value = FALSE)),
                    column(4, strong(nm)),
                    column(2, span(df$size_fmt[i],  class = "size-col")),
                    column(2, span(df$class[i],     class = "class-col"))
                )
            })

            tagList(header, rows)
        })

        observeEvent(input$remove_btn, {
            df <- obj_data()
            selected <- df$name[vapply(df$name, function(nm) isTRUE(input[[paste0("chk_", nm)]]), logical(1))]

            if (length(selected) == 0) {
                showNotification("No objects selected.", type = "warning")
                return()
            }

            rm(list = selected, envir = env)
            gc()
            obj_data(get_obj_table())
            showNotification(
                paste0("Removed ", length(selected), " object(s): ", paste(selected, collapse = ", ")),
                type = "message", duration = 5
            )
        })

        observeEvent(input$select_all, {
            df <- obj_data()
            for (nm in df$name) updateCheckboxInput(session, paste0("chk_", nm), value = TRUE)
        })

        observeEvent(input$desel_all, {
            df <- obj_data()
            for (nm in df$name) updateCheckboxInput(session, paste0("chk_", nm), value = FALSE)
        })
    }

    shinyApp(ui, server)
}
