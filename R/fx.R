# fx
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 2632-2661). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

fx <- function(df, title = "your_data") {
  
  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("fx() needs openxlsx for Excel output. install.packages('openxlsx').")
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
  addWorksheet(wb, "Sheet 1")
  writeData(wb, "Sheet 1", df)
  setColWidths(wb, "Sheet 1", cols = 1:ncol(df), widths = "auto")
  freezePane(wb, "Sheet 1", firstRow = TRUE)
  addFilter(wb, "Sheet 1", rows = 1, cols = 1:ncol(df))


  saveWorkbook(wb, excel_filename, overwrite = TRUE)

}
