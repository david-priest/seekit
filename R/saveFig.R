# saveFig.R
#
# Replacement for f2().  Saves a plot object (ggplot or base_plot) to PDF,
# and optionally saves the underlying plot data to Excel and the plot object
# itself to RData.  Also keeps a rolling copy of the source RMD and
# sessionInfo in the output folder.
#
# USAGE
#   source("R/saveFig.R")
#   p <- ggplot(...)          # ggplot object  — OR —
#   p <- plotAbundancesBase() # base_plot object
#   saveFig(p, w = 20, h = 10, title = "My figure")
#   saveFig(p, w = 20, h = 10, title = "My figure", saveExcel = TRUE)
#   saveFig(p, w = 20, h = 10, title = "My figure", format = "svg")
#   saveFig(p, w = 20, h = 10, title = "My figure", format = "svg",
#           svg_mode = "grouped_auto")

.countTag <- function(pattern, x) {
  m <- gregexpr(pattern, x, perl = TRUE)[[1L]]
  if (length(m) == 1L && m[1L] == -1L) 0L else length(m)
}

.findGroupEnd <- function(lines, start_idx) {
  depth <- 0L
  n <- length(lines)
  for (i in start_idx:n) {
    depth <- depth + .countTag("<g\\b", lines[i]) - .countTag("</g>", lines[i])
    if (depth <= 0L) return(i)
  }
  NA_integer_
}

.extractClipIdFromG <- function(line) {
  m <- regexec("<g[^>]*clip-path=['\"]url\\(#([^)'\"]+)\\)['\"][^>]*>", line, perl = TRUE)
  r <- regmatches(line, m)[[1L]]
  if (length(r) >= 2L) r[2L] else NA_character_
}

.discoverClipKinds <- function(lines) {
  panel_ids  <- character()
  canvas_ids <- character()

  n <- length(lines)
  if (n < 3L) return(list(panel_ids = panel_ids, canvas_ids = canvas_ids))

  id_pat <- "^\\s*<clipPath id=['\"]([^'\"]+)['\"]>\\s*$"
  rect_pat <- paste0(
    "^\\s*<rect\\s+x=['\"]([0-9.+-eE]+)['\"]\\s+y=['\"]([0-9.+-eE]+)['\"]",
    "\\s+width=['\"]([0-9.+-eE]+)['\"]\\s+height=['\"]([0-9.+-eE]+)['\"][^>]*>\\s*$"
  )

  for (i in seq_len(n - 2L)) {
    id_m <- regexec(id_pat, lines[i], perl = TRUE)
    id_r <- regmatches(lines[i], id_m)[[1L]]
    if (length(id_r) < 2L) next

    rect_m <- regexec(rect_pat, lines[i + 1L], perl = TRUE)
    rect_r <- regmatches(lines[i + 1L], rect_m)[[1L]]
    if (length(rect_r) < 5L) next

    id <- id_r[2L]
    x  <- suppressWarnings(as.numeric(rect_r[2L]))
    y  <- suppressWarnings(as.numeric(rect_r[3L]))

    if (!is.finite(x) || !is.finite(y)) next
    if (abs(x) < 1e-8 && abs(y) < 1e-8) {
      canvas_ids <- c(canvas_ids, id)
    } else {
      panel_ids <- c(panel_ids, id)
    }
  }

  list(panel_ids = unique(panel_ids), canvas_ids = unique(canvas_ids))
}

.groupSvgPanelsAuto <- function(lines, panel_prefix = "panel") {
  clips <- .discoverClipKinds(lines)
  panel_ids <- clips$panel_ids
  canvas_ids <- clips$canvas_ids

  if (length(panel_ids) == 0L) return(lines)

  out <- character()
  i <- 1L
  n <- length(lines)
  panel_idx <- 0L

  while (i <= n) {
    line <- lines[i]
    clip_id <- .extractClipIdFromG(line)

    if (!is.na(clip_id) && clip_id %in% panel_ids) {
      j <- .findGroupEnd(lines, i)
      if (is.na(j)) {
        out <- c(out, line)
        i <- i + 1L
        next
      }

      # In svglite output, axis/text for a panel are usually in the next
      # sibling <g> clipped to the full canvas. Keep both together.
      l <- j
      k <- j + 1L
      while (k <= n && grepl("^\\s*$", lines[k])) k <- k + 1L
      if (k <= n) {
        clip_id2 <- .extractClipIdFromG(lines[k])
        if (!is.na(clip_id2) && clip_id2 %in% canvas_ids) {
          l2 <- .findGroupEnd(lines, k)
          if (!is.na(l2)) l <- l2
        }
      }

      panel_idx <- panel_idx + 1L
      out <- c(
        out,
        sprintf('<g id="%s_%03d">', panel_prefix, panel_idx),
        lines[i:l],
        "</g>"
      )
      i <- l + 1L
      next
    }

    out <- c(out, line)
    i <- i + 1L
  }

  out
}

# Quick annotation:
# - format = "pdf": standard PDF export.
# - format = "svg", svg_mode = "plain": fast ungrouped SVG from svglite.
# - format = "svg", svg_mode = "grouped_auto": post-processes SVG and wraps
#   detected panel blocks in <g id="panel_###"> groups for easier editing.
saveFig <- function(p,
                    w         = 20,
                    h         = 20,
                    title     = "your_title",
                    saveRData = FALSE,
                    saveExcel = FALSE,
                    output_root = NULL,
                    format    = c("pdf", "svg"),
                    svg_mode  = c("plain", "grouped_auto")) {
  # format:  output file format
  # svg_mode: for format='svg', choose plain (fast) or grouped_auto (panel groups)

  format   <- match.arg(format)
  svg_mode <- match.arg(svg_mode)

  # ── Output folder ────────────────────────────────────────────────────────────
  # Source qmd/Rmd name without a hard rmdhelp dependency (rmdhelp is optional).
  rmd_file_path      <- .wl_source_file()
  rmd_file_name      <- if (is.na(rmd_file_path)) "figures" else basename(rmd_file_path)
  rmd_file_name_noext <- tools::file_path_sans_ext(rmd_file_name)

  folder_base <- if (is.null(output_root))
                   paste0(Sys.Date(), "_", rmd_file_name_noext)
                 else
                   file.path(output_root, paste0(Sys.Date(), "_", rmd_file_name_noext))
  folder_name <- folder_base
  if (!dir.exists(folder_name)) dir.create(folder_name, recursive = TRUE)

  timestamp <- format(Sys.time(), "%Y-%m-%d_%H.%M.%S")

  # ── Plot export ─────────────────────────────────────────────────────────────
  fig_filename <- file.path(folder_name, paste0(timestamp, "_", title, ".", format))

  if (format == "pdf") {
    grDevices::pdf(fig_filename, width = w, height = h)
    invisible(print(p))
    dev.off()
  } else {
    if (!requireNamespace("svglite", quietly = TRUE)) {
      stop("Package 'svglite' is required for SVG export. Please install it.")
    }

    if (svg_mode == "plain") {
      svglite::svglite(fig_filename, width = w, height = h)
      invisible(print(p))
      dev.off()
    } else {
      tmp_svg <- tempfile(pattern = "savefig_plain_", fileext = ".svg")
      svglite::svglite(tmp_svg, width = w, height = h)
      invisible(print(p))
      dev.off()

      svg_lines <- readLines(tmp_svg, warn = FALSE)
      svg_lines <- .groupSvgPanelsAuto(svg_lines, panel_prefix = "panel")
      writeLines(svg_lines, fig_filename, useBytes = TRUE)
      unlink(tmp_svg)
    }
  }

  # ── RData ────────────────────────────────────────────────────────────────────
  if (saveRData) {
    rdata_filename <- file.path(folder_name, paste0(timestamp, "_", title, ".RData"))
    cat("Saving RData...\n")
    save(p, file = rdata_filename)
  }

  # ── Excel (plot data) ────────────────────────────────────────────────────────
  if (saveExcel && !requireNamespace("openxlsx", quietly = TRUE)) {
    warning("openxlsx is not installed; skipping Excel data save. install.packages('openxlsx').")
  } else if (saveExcel) {
    library(openxlsx)

    # Extract the underlying data regardless of plot type:
    #   ggplot  → p$data
    #   base_plot (plotAbundancesBase / plotAbundancesSimple) → p$data
    plot_data <- NULL
    if (!is.null(p$data) && is.data.frame(p$data)) {
      plot_data <- p$data
    }

    if (is.null(plot_data)) {
      warning("saveExcel = TRUE but no data.frame found in p$data — Excel file not written.")
    } else {
      excel_filename <- file.path(folder_name,
                                  paste0(timestamp, "_", title, "_data.xlsx"))
      cat("Saving Excel...\n")

      wb <- createWorkbook()
      modifyBaseFont(wb, fontSize = 12, fontColour = "black", fontName = "Arial")
      addWorksheet(wb, "Plot data")
      writeData(wb, "Plot data", plot_data)
      setColWidths(wb, "Plot data", cols = seq_len(ncol(plot_data)), widths = "auto")
      freezePane(wb,  "Plot data", firstRow = TRUE)
      addFilter(wb,   "Plot data", rows = 1, cols = seq_len(ncol(plot_data)))
      saveWorkbook(wb, excel_filename, overwrite = TRUE)
    }
  }

  # ── Rolling RMD + sessionInfo copy (no more than once per 20 minutes) ────────
  recent_files    <- list.files(folder_name,
                                pattern    = "\\.rmd$|_sessioninfo\\.txt$",
                                full.names = TRUE,
                                ignore.case = TRUE)
  latest_mtime    <- if (length(recent_files) > 0)
                       max(file.info(recent_files)$mtime, na.rm = TRUE)
                     else NA
  needs_rmd_copy  <- is.na(latest_mtime) ||
                     difftime(Sys.time(), latest_mtime, units = "mins") > 20

  if (needs_rmd_copy) {
    if (!is.na(rmd_file_path) && file.exists(rmd_file_path))
      file.copy(from = rmd_file_path,
                to   = file.path(folder_name,
                                 paste0(timestamp, "_", rmd_file_name)))

    writeLines(capture.output(sessionInfo()),
               con = file.path(folder_name,
                               paste0(timestamp, "_sessioninfo.txt")))

    cat("Figure, RMD file, and session information saved in folder",
        folder_name, "\n")
  } else {
    cat("Figure saved in folder", folder_name,
        "— RMD/sessionInfo skipped (recent copy exists).\n")
  }

  invisible(NULL)
}
