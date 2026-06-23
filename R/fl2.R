# fl2
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 2547-2626). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

fl2 <- function(p_list, w = 20, h = 20, title = "your_title", saveRData = FALSE, saveExcel = FALSE) {
  # Name of the source qmd/Rmd, resolved without a hard rmdhelp dependency
  # (rmdhelp is GitHub-only and optional). Falls back to "figures" if unknown.
  rmd_file_path <- .wl_source_file()
  rmd_file_name <- if (is.na(rmd_file_path)) "figures" else basename(rmd_file_path)
  rmd_file_name_noext <- tools::file_path_sans_ext(rmd_file_name)

  # Anchor folder to here::here() (project root) rather than getwd().
  # See f2.R for full rationale.
  folder_name <- here::here(paste0(Sys.Date(), "_", rmd_file_name_noext))

  # Check if the directory exists, and create it if it doesn't
  if (!dir.exists(folder_name)) {
    dir.create(folder_name, recursive = TRUE)
  }
  
  # Create the timestamp for the file names
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H.%M.%S")
  
  # Create the filename for the PDF figure
  pdf_filename <- paste0(folder_name, "/", timestamp, "_", title, ".pdf")
  
  # Cairo banned: always base grDevices::pdf() (matches f()/f2()). Trade-off: some
  # Unicode glyphs (e.g. the "∗" math asterisk) render as fallback dots — use ASCII
  # "*" in labels, as elsewhere in the codebase.
  grDevices::pdf(pdf_filename, width = w, height = h)

  # Loop through each figure in p_list and print it to a new page in the PDF
  for(p in p_list) {
    invisible(print(p))
  }

  dev.off()

  # Save RData
  if (saveRData) {
    rdata_filename <- file.path(folder_name, paste0(timestamp, "_", title, ".RData"))
    cat("Saving RData...\n")
    save(p_list, file = rdata_filename)
  }

  # Save Excel
  if (saveExcel && !requireNamespace("openxlsx", quietly = TRUE)) {
    warning("openxlsx is not installed; skipping Excel data save. install.packages('openxlsx').")
  } else if (saveExcel) {
    library(openxlsx)
    excel_filename <- paste0(folder_name, "/", timestamp, "_", title, "_data.xlsx")
    cat("Saving Excel...\n")
    wb <- createWorkbook()
    modifyBaseFont(wb, fontSize = 12, fontColour = "black", fontName = "Arial")
    
    for (i in seq_along(p_list)) {
      p <- p_list[[i]]
      sheet_name <- names(p_list)[i]
      if (!is.null(p$data)) {
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet_name, p$data)
        setColWidths(wb, sheet_name, cols = 1:ncol(p$data), widths = "auto")
        freezePane(wb, sheet_name, firstRow = TRUE)
        addFilter(wb, sheet_name, rows = 1, cols = 1:ncol(p$data))
      }
    }
    
    saveWorkbook(wb, excel_filename, overwrite = TRUE)
  }

  # Check if an RMD or sessionInfo file has been saved within the last 5 minutes
  recent_files <- list.files(folder_name, pattern = "\\.rmd$|_sessioninfo\\.txt$", full.names = TRUE)
  latest_file_time <- max(file.info(recent_files)$mtime, na.rm = TRUE)
  current_time <- Sys.time()

  # If no RMD or sessionInfo file saved in the last 5 minutes, save a copy
  if (is.na(latest_file_time) || difftime(current_time, latest_file_time, units = "mins") > 20) {
    if (!is.na(rmd_file_path) && file.exists(rmd_file_path)) {
      rmd_filename <- paste0(folder_name, "/", timestamp, "_", rmd_file_name)
      file.copy(from = rmd_file_path, to = rmd_filename)
    }

    # Save the output of sessionInfo() as a text file
    sessioninfo_filename <- paste0(folder_name, "/", timestamp, "_sessioninfo.txt")
    writeLines(capture.output(sessionInfo()), con = sessioninfo_filename)

    cat("Figures, RMD file, and session information saved successfully in folder", folder_name, "\n")
  } else {
    cat("Figures saved successfully in folder", folder_name, ". RMD file and session information were not saved as a recent copy already exists.\n")
  }
}
