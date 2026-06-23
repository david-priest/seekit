# f2
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 2418-2484). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

f2 <- function(p, w = 20, h = 20, title = "your_title",
               saveRData = FALSE, saveExcel = FALSE,
               format = c("pdf", "png"), png_res = 300) {
  # Format selection. "pdf" stays the default; "png" produces a hi-res raster
  # (ragg::agg_png at `png_res` DPI if ragg is installed, falling back to
  # grDevices::png with type = "cairo"). w/h continue to be interpreted as
  # inches for both formats.
  format <- match.arg(format)
  if (format == "png" && png_res <= 0) {
    stop("png_res must be a positive integer (DPI). Got: ", png_res)
  }

  # Get the name of the current RMD file
  library(rmdhelp)
  rmd_file_path <- get_this_rmd_file()
  rmd_file_name <- basename(rmd_file_path)
  rmd_file_name_noext <- tools::file_path_sans_ext(rmd_file_name)

  # Anchor the dated output folder to the project root (the directory
  # containing the .here marker) via here::here(), NOT to getwd(). In
  # Positron / RStudio shared sessions getwd() can drift to another
  # project's dir even after setwd() runs in setup_libraries -- the
  # original relative-path code would then write the figure into that
  # other project's directory. Using here::here() makes the folder
  # location robust to that drift. (G13 trap from the CMV plotting ledger.)
  folder_name <- here::here(paste0(Sys.Date(), "_", rmd_file_name_noext))

  # Check if the directory exists, and create it if it doesn't
  if (!dir.exists(folder_name)) {
    dir.create(folder_name, recursive = TRUE)
  }

  # Create the timestamp for the file names
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H.%M.%S")

  # Create the filename for the figure (extension follows `format`).
  fig_filename <- file.path(folder_name,
                            paste0(timestamp, "_", title, ".", format))

  if (format == "pdf") {
    # Use base grDevices::pdf() rather than cairo_pdf.
    # Why: cairo_pdf applies OpenType ligature features via Pango, so titles
    # containing "fi" / "fl" etc. (e.g. "Non-specific") get rendered with the
    # ligature glyph (ﬁ). Illustrator can't decompose that glyph back into
    # individual characters within a single text run -> the title gets split
    # into multiple unrelated text objects ("Non-speci" + glyph + "c"), which
    # breaks editability. Base pdf() doesn't apply ligatures so titles stay as
    # single editable strings in Illustrator. Trade-off: some Unicode glyphs
    # (e.g. U+2217 math asterisk) may render as fallback dots if the Type 1
    # font doesn't cover them; we accept this since the codebase's significance
    # asterisks are plain ASCII `*` (U+002A), not the math asterisk.
    grDevices::pdf(fig_filename, width = w, height = h)
  } else {
    # Hi-res PNG via ragg::agg_png (Cairo-free, crisp text antialiasing). We do
    # NOT fall back to grDevices::png(type = "cairo") — Cairo is banned (it fails
    # / segfaults on large plots). w/h in inches, scaled by png_res to pixels.
    if (!requireNamespace("ragg", quietly = TRUE))
      stop("ragg is required for Cairo-free PNG output. install.packages('ragg').")
    ragg::agg_png(fig_filename, width = w, height = h,
                  units = "in", res = png_res)
  }
  invisible(print(p))
  dev.off()

  # ── Provenance: write a MANIFEST.jsonl entry per save ──────────────────────
  # One JSON object per line, appended. Captures git commit, source qmd,
  # chunk label, timestamp, and figure path so any saved figure can be
  # traced back to the code that produced it (figures get .gitignore'd
  # but MANIFEST stays in the repo). Backward compatible: degrades
  # silently if gert / jsonlite / knitr aren't installed.
  manifest_path <- file.path(folder_name, "MANIFEST.jsonl")
  prov <- list(
    fig         = basename(fig_filename),
    title       = title,
    fig_format  = format,
    width_in    = w,
    height_in   = h,
    timestamp   = timestamp,
    saved_at    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    qmd_path    = tryCatch(rmd_file_path, error = function(e) NA_character_),
    chunk_label = tryCatch(knitr::opts_current$get("label"),
                           error = function(e) NA_character_),
    git_commit  = tryCatch(
      if (requireNamespace("gert", quietly = TRUE) &&
          !inherits(gi <- gert::git_info(repo = dirname(rmd_file_path)),
                    "try-error")) substr(gi$commit, 1, 12) else NA_character_,
      error = function(e) NA_character_
    ),
    git_branch  = tryCatch(
      if (requireNamespace("gert", quietly = TRUE) &&
          !inherits(gi2 <- gert::git_info(repo = dirname(rmd_file_path)),
                    "try-error")) gi2$shorthand else NA_character_,
      error = function(e) NA_character_
    ),
    r_version   = paste(R.version$major, R.version$minor, sep = ".")
  )
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    tryCatch({
      jline <- jsonlite::toJSON(prov, auto_unbox = TRUE, null = "null",
                                na = "null")
      cat(as.character(jline), "\n", file = manifest_path, append = TRUE,
          sep = "")
    }, error = function(e) invisible(NULL))
  }

 # Save RData
  if (saveRData) {
    rdata_filename <- file.path(folder_name, paste0(timestamp, "_", title, ".RData"))
    cat("Saving RData...\n")
    save(p, file = rdata_filename)
  }

  # Save Excel
  #
  # Robust to patchwork composites: if p is a patchwork, we recurse through the
  # component plots, write each plot's main $data to its own sheet, and *also*
  # scan each plot's layers for stat-bracket data (group1/group2/label/y.position
  # style columns) which is what plotAbundancesDiffStrip's stats strip carries.
  # This way you get the data table AND the stats table in one xlsx, even when
  # the figure is `strip / data` glued together with patchwork.
  if (saveExcel) {
    library(openxlsx)
    excel_filename <- file.path(folder_name, paste0(timestamp, "_", title, "_data.xlsx"))
    cat("Saving Excel...\n")

    # --- helpers ---------------------------------------------------------------
    # Recursively flatten a patchwork into a list of ggplot objects, in
    # roughly the order they were combined.
    collect_plots <- function(x) {
      if (inherits(x, "patchwork")) {
        out <- list()
        subs <- tryCatch(x$patches$plots, error = function(e) NULL)
        if (length(subs)) for (sub in subs) out <- c(out, collect_plots(sub))
        # The patchwork object itself also carries the "last" ggplot's slots
        # (data, layers, mapping, ...). Strip the patchwork class to access
        # them as a plain ggplot.
        main <- x
        class(main) <- setdiff(class(main), "patchwork")
        out <- c(out, list(main))
        return(out)
      }
      if (inherits(x, "ggplot")) return(list(x))
      list()
    }

    # Pull a plot's primary data, falling back to first non-waiver layer.
    get_plot_data <- function(plt) {
      d <- plt$data
      if (!inherits(d, "waiver") && !is.null(d) && length(d) > 0) {
        df <- try(as.data.frame(d), silent = TRUE)
        if (!inherits(df, "try-error") && nrow(df) > 0) return(df)
      }
      for (l in plt$layers) {
        ld <- l$data
        if (!inherits(ld, "waiver") && !is.null(ld)) {
          df <- try(as.data.frame(ld), silent = TRUE)
          if (!inherits(df, "try-error") && nrow(df) > 0) return(df)
        }
      }
      NULL
    }

    # Heuristic: a stats layer carries one of these columns.
    stats_cols <- c("p.adj.signif", "p.signif", "p.adj", "p.value",
                    "group1", "group2", "y.position", "label", "xmin", "xmax")
    get_stats_data <- function(plt, main_df) {
      out <- list()
      for (l in plt$layers) {
        ld <- l$data
        if (inherits(ld, "waiver") || is.null(ld)) next
        df <- try(as.data.frame(ld), silent = TRUE)
        if (inherits(df, "try-error") || nrow(df) == 0) next
        # Skip if it's identical to the plot's main data
        if (!is.null(main_df) &&
            identical(dim(df), dim(main_df)) &&
            identical(sort(names(df)), sort(names(main_df)))) next
        if (any(stats_cols %in% names(df))) out[[length(out) + 1L]] <- df
      }
      out
    }

    sanitize_sheet <- function(s) {
      s <- gsub("[\\\\/?*\\[\\]:]", "_", s, perl = TRUE)
      substr(s, 1, 31)
    }

    # --- collect plots & write -------------------------------------------------
    plots <- collect_plots(p)
    if (length(plots) == 0) {
      cat("No ggplot/patchwork content to save; skipping Excel.\n")
    } else {
      wb <- createWorkbook()
      modifyBaseFont(wb, fontSize = 12, fontColour = "black", fontName = "Arial")
      any_written <- FALSE

      add_sheet <- function(name, df) {
        sn <- sanitize_sheet(name)
        # Ensure uniqueness against existing sheet names
        existing <- names(wb)
        if (sn %in% existing) {
          k <- 2L
          while (paste0(substr(sn, 1, 28), "_", k) %in% existing) k <- k + 1L
          sn <- paste0(substr(sn, 1, 28), "_", k)
        }
        addWorksheet(wb, sn)
        writeData(wb, sn, df)
        if (ncol(df) >= 1) setColWidths(wb, sn, cols = 1:ncol(df), widths = "auto")
        freezePane(wb, sn, firstRow = TRUE)
        if (ncol(df) >= 1) addFilter(wb, sn, rows = 1, cols = 1:ncol(df))
      }

      # --- Preferred path: long-format data + stats stashed as attributes ---
      # plotAbundancesDiffStrip() and friends attach `source_data`,
      # `source_stats`, and `n_per_condition` to the returned patchwork.
      # If present, write them as their own sheets and skip the per-plot
      # extraction (which would otherwise spread data across N sheets and
      # lose stats info because the rendered stats layers only carry x/y/
      # label/vj).
      src_data  <- attr(p, "source_data")
      src_stats <- attr(p, "source_stats")
      src_n     <- attr(p, "n_per_condition")
      used_attrs <- FALSE
      if (!is.null(src_data)) {
        df_src <- try(as.data.frame(src_data), silent = TRUE)
        if (!inherits(df_src, "try-error") && nrow(df_src) > 0) {
          add_sheet("data", df_src)
          any_written <- TRUE
          used_attrs  <- TRUE
        }
      }
      if (!is.null(src_stats)) {
        df_stats <- try(as.data.frame(src_stats), silent = TRUE)
        if (!inherits(df_stats, "try-error") && nrow(df_stats) > 0) {
          add_sheet("stats", df_stats)
          any_written <- TRUE
          used_attrs  <- TRUE
        }
      }
      if (!is.null(src_n)) {
        df_n <- try(as.data.frame(src_n), silent = TRUE)
        if (!inherits(df_n, "try-error") && nrow(df_n) > 0) {
          add_sheet("n_per_condition", df_n)
          any_written <- TRUE
          used_attrs  <- TRUE
        }
      }

      # --- Fallback path: per-plot extraction (legacy / generic ggplots) ---
      if (!used_attrs) {
        # Single-plot case (no patchwork) -> keep legacy "Sheet 1" name.
        single_plot <- length(plots) == 1L
        for (i in seq_along(plots)) {
          plt <- plots[[i]]
          main_df <- get_plot_data(plt)
          prefix <- if (single_plot) "Sheet 1" else paste0("Plot_", i)
          if (!is.null(main_df)) {
            add_sheet(prefix, main_df)
            any_written <- TRUE
          }
          stats_list <- get_stats_data(plt, main_df)
          for (j in seq_along(stats_list)) {
            nm <- if (length(stats_list) == 1L) paste0(prefix, "_stats")
                  else paste0(prefix, "_stats_", j)
            add_sheet(nm, stats_list[[j]])
            any_written <- TRUE
          }
        }
      }

      if (any_written) {
        saveWorkbook(wb, excel_filename, overwrite = TRUE)
      } else {
        cat("No tabular data found in plot; skipping Excel.\n")
      }
    }
  }

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
