# fxl
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 2667-2712). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

fxl <- function(df_list, title = "your_data", separate_sheets = TRUE) {
  
  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("fxl() needs openxlsx for Excel output. install.packages('openxlsx').")
  library(openxlsx)

  # Source qmd/Rmd name without a hard rmdhelp dependency (rmdhelp is optional).
  rmd_file_path <- .wl_source_file()
  rmd_file_name <- if (is.na(rmd_file_path)) "figures" else basename(rmd_file_path)
  rmd_file_name_noext <- tools::file_path_sans_ext(rmd_file_name)

  # Anchor folder to here::here() (project root) rather than getwd().
  # See f2.R for full rationale.
  folder_name <- here::here(paste0(Sys.Date(), "_", rmd_file_name_noext))

  if (!dir.exists(folder_name)) {
    dir.create(folder_name, recursive = TRUE)
  }

  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  excel_filename <- file.path(folder_name, paste0(timestamp, "_", title, ".xlsx"))
  
  wb <- createWorkbook()
  modifyBaseFont(wb, fontSize = 12, fontColour = "black", fontName = "Arial")
  
  if (separate_sheets) {
    # Use the list element names as sheet names when the list is named.
    # Falls back to "Sheet i" for any unnamed / blank entry. Sanitises for
    # Excel's sheet-name rules: no `\ / ? * [ ] :`, max 31 chars; appends a
    # suffix on collision so we don't error mid-write.
    sanitize_sheet <- function(s) {
      s <- gsub("[\\\\/?*\\[\\]:]", "_", s, perl = TRUE)
      substr(s, 1, 31)
    }
    list_names <- names(df_list)
    if (is.null(list_names)) list_names <- rep("", length(df_list))

    for (i in seq_along(df_list)) {
      df <- df_list[[i]]
      raw_name  <- list_names[i]
      base_name <- if (!is.na(raw_name) && nzchar(raw_name)) raw_name else paste0("Sheet ", i)
      sheet_name <- sanitize_sheet(base_name)

      existing <- names(wb)
      if (sheet_name %in% existing) {
        k <- 2L
        while (paste0(substr(sheet_name, 1, 28), "_", k) %in% existing) k <- k + 1L
        sheet_name <- paste0(substr(sheet_name, 1, 28), "_", k)
      }

      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, df)
      if (ncol(df) >= 1) setColWidths(wb, sheet_name, cols = 1:ncol(df), widths = "auto")
      freezePane(wb, sheet_name, firstRow = TRUE)
      if (ncol(df) >= 1) addFilter(wb, sheet_name, rows = 1, cols = 1:ncol(df))
    }
  } else {
    addWorksheet(wb, "Concatenated Data")
    start_row <- 1
    for (i in seq_along(df_list)) {
      df <- df_list[[i]]
      writeData(wb, "Concatenated Data", df, startRow = start_row)
      setColWidths(wb, "Concatenated Data", cols = 1:ncol(df), widths = "auto")
      freezePane(wb, "Concatenated Data", firstRow = TRUE)
      addFilter(wb, "Concatenated Data", rows = start_row, cols = 1:ncol(df))
      start_row <- start_row + nrow(df) + 2  # Move to the next row after the current data frame plus one empty row
    }
  }

  saveWorkbook(wb, excel_filename, overwrite = TRUE)
  cat("Data frames saved successfully in", excel_filename, "\n")
}
