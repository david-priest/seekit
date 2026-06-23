# f
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 2368-2412). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

f <- function(p, w = 20, h = 20, title = "your_title", cairo = FALSE) {
  # Get the name of the current RMD file
  library(rmdhelp)
  rmd_file_path <- get_this_rmd_file()
  rmd_file_name <- basename(rmd_file_path)
  rmd_file_name_noext <- tools::file_path_sans_ext(rmd_file_name)

  # Anchor the dated output folder to here::here() (project root) rather
  # than to getwd(). See f2.R for full rationale.
  folder_name <- here::here(paste0(Sys.Date(), "_", rmd_file_name_noext))

  # Check if the directory exists, and create it if it doesn't
  if (!dir.exists(folder_name)) {
    dir.create(folder_name, recursive = TRUE)
  }
  
  # Create the timestamp for the file names
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H.%M.%S")
  
  # Create the filename for the PDF figure using file.path for compatibility
  pdf_filename <- file.path(folder_name, paste0(timestamp, "_", title, ".pdf"))
  
  # Default to base grDevices::pdf() so every text element (facet/strip titles,
  # axis labels, legend text) is emitted as a SEPARATE PDF text object —
  # Illustrator then imports each as its own editable text box. cairo_pdf()
  # groups adjacent text runs into one text frame (and applies OpenType
  # ligatures, splitting "fi"/"fl" titles), so multiple panel titles arrive
  # merged into a single, hard-to-edit box. See the same rationale in f2.R.
  #
  # Trade-off: base pdf() can't render non-Latin-1 glyphs such as the U+2217
  # math asterisk "∗" (renders as a fallback dot). Pass `cairo = TRUE` only when
  # you need such Unicode glyphs and aren't going to edit the text in
  # Illustrator. Significance marks in this codebase use the ASCII "*" (U+002A),
  # which base pdf() renders fine.
  # Cairo is banned (cairo_pdf fails/segfaults + applies ligatures that break
  # Illustrator editability). Always base grDevices::pdf(); the `cairo` arg is
  # retained for back-compat but ignored.
  grDevices::pdf(pdf_filename, width = w, height = h)
  invisible(print(p))
  dev.off()

  # Check if an RMD or sessionInfo file has been saved within the last 5 minutes
  recent_files <- list.files(folder_name, pattern = "\\.rmd$|_sessioninfo\\.txt$", full.names = TRUE)
  latest_file_time <- max(file.info(recent_files)$mtime, na.rm = TRUE)
  current_time <- Sys.time()

  # If no RMD or sessionInfo file saved in the last 5 minutes, save a copy
  if (is.na(latest_file_time) || difftime(current_time, latest_file_time, units = "mins") > 20) {
    rmd_filename <- file.path(folder_name, paste0(timestamp, "_", rmd_file_name))
    file.copy(from = rmd_file_path, to = rmd_filename)

    # Save the output of sessionInfo() as a text file
    sessioninfo_filename <- file.path(folder_name, paste0(timestamp, "_sessioninfo.txt"))
    writeLines(capture.output(sessionInfo()), con = sessioninfo_filename)

    cat("Figure, RMD file, and session information saved successfully in folder", folder_name, "\n")
  } else {
    cat("Figure saved successfully in folder", folder_name, ". RMD file and session information were not saved as a recent copy already exists.\n")
  }
}
