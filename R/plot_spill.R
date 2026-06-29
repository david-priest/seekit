#' @title Plot per-interaction spillover scatter panels
#' @description Build a grid of biaxial scatter panels for the channel pairs in a
#'   spillover matrix, to eyeball compensation. Supports two modes: `"spill_matrix"`
#'   (every non-zero entry) and `"valid_interactions"` (only physically plausible
#'   interactions: same element, +16 oxide, +/-1 isotope), and groups/orders panels
#'   by source or destination channel. Panels are titled with the spillover
#'   coefficient and the source/dest antigen names from `panel_table`.
#' @param spill_matrix Numeric spillover matrix, channel names on both dims
#'   (rows = source, cols = destination).
#' @param cytof_data A `SingleCellExperiment` (CATALYST-style) with an `exprs`
#'   assay and a `sample_id` colData column; antigens as rownames.
#' @param panel_table Data frame mapping channels to antigens via `fcs_colname`
#'   and `antigen` columns.
#' @param n_cells Cells to sub-sample per sample for the scatters, Default: 1000.
#' @param my_seed RNG seed for the sub-sample, Default: 1234.
#' @param exclude_channels Optional character vector of source channels to skip.
#' @param n_col Columns in the final panel grid, Default: 6.
#' @param mode One of `"spill_matrix"` / `"valid_interactions"`, Default: "spill_matrix".
#' @param valid_interactions Which interaction types count as valid in
#'   `"valid_interactions"` mode, Default: c("same","ox","one_plus","one_minus").
#' @param rasterize,dpi Retained for back-compatibility; rasterization is a no-op
#'   here (the package is Cairo-/ggrastr-free — hexbin panels render fine as vector).
#' @param title_size,axis_size Panel title / axis text sizes, Default: 8 / 7.
#' @param organize_by Group + order panels by `"source"` or `"destination"`,
#'   Default: "source".
#' @return A composed grid (`cowplot::plot_grid`) of scatter panels, or `NULL` if
#'   there is nothing to plot.
#' @note Relies on a `plotScatter()` (the CATALYST biaxial plotter) being available
#'   in scope, as the analysis workflow attaches it.
#' @export
#' @importFrom cowplot plot_grid
#' @importFrom readr parse_number
plot_spill <- function(spill_matrix, cytof_data, panel_table, n_cells = 1000,
                       my_seed = 1234, exclude_channels = NULL, n_col = 6,
                       mode = c("spill_matrix", "valid_interactions"),
                       valid_interactions = c("same", "ox", "one_plus", "one_minus"),
                       rasterize = TRUE, dpi = 150, title_size = 8, axis_size = 7,
                       organize_by = c("source", "destination")) {
  mode <- match.arg(mode)
  organize_by <- match.arg(organize_by)

  set.seed(my_seed)
  idx <- split(seq(ncol(cytof_data)), cytof_data$sample_id)
  idx <- lapply(idx, function(.) sample(., min(n_cells, length(.))))
  cytof_data <- cytof_data[, unlist(idx)]

  is_valid_interaction <- function(source, dest, allowed = valid_interactions) {
    same_element <- (gsub("[0-9]", "", dest) == gsub("[0-9]", "", source))
    ox         <- (readr::parse_number(dest) == readr::parse_number(source) + 16)
    one_plus   <- (readr::parse_number(dest) == readr::parse_number(source) + 1)
    one_minus  <- (readr::parse_number(dest) == readr::parse_number(source) - 1)
    (("same" %in% allowed && same_element) || ("ox" %in% allowed && ox) ||
       ("one_plus" %in% allowed && one_plus) || ("one_minus" %in% allowed && one_minus))
  }

  if (mode == "spill_matrix") {
    relevant_entries <- which(spill_matrix != 0, arr.ind = TRUE)
    if (nrow(relevant_entries) == 0) {
      warning("No non-zero entries in the spill matrix. No plots to generate.")
      return(NULL)
    }
  } else {
    relevant_entries <- which(spill_matrix >= 0 & spill_matrix <= 1, arr.ind = TRUE)
    if (nrow(relevant_entries) == 0) {
      warning("No valid entries in the spill matrix. No plots to generate.")
      return(NULL)
    }
  }

  if (organize_by == "source") {
    order_channels <- rownames(spill_matrix)[relevant_entries[, "row"]]
  } else {
    order_channels <- colnames(spill_matrix)[relevant_entries[, "col"]]
  }
  order_numbers <- as.numeric(gsub("[^0-9]", "", order_channels))
  relevant_entries <- relevant_entries[order(order_numbers), , drop = FALSE]

  grouped_plots <- list()
  for (i in seq_len(nrow(relevant_entries))) {
    source <- rownames(spill_matrix)[relevant_entries[i, "row"]]
    dest   <- colnames(spill_matrix)[relevant_entries[i, "col"]]
    if (source == dest) next
    if (!is.null(exclude_channels) && source %in% exclude_channels) next

    if (mode == "valid_interactions") {
      if (!is_valid_interaction(source, dest)) next
      if ("nonzero" %in% valid_interactions && spill_matrix[source, dest] == 0) next
    }

    source_antigen <- panel_table$antigen[panel_table$fcs_colname == source]
    dest_antigen   <- panel_table$antigen[panel_table$fcs_colname == dest]
    if (length(source_antigen) < 1 || length(dest_antigen) < 1) {
      message(paste0("Could not map channels ", source, " and/or ", dest,
                     " to antigens. Skipping these channels."))
      next
    }
    if (!(source_antigen %in% rownames(cytof_data)) || !(dest_antigen %in% rownames(cytof_data))) {
      message(paste0("Antigens ", source_antigen, " and/or ", dest_antigen,
                     " are not valid in the CyTOF data. Skipping these antigens."))
      next
    }

    spill_value <- sprintf("%.4f", spill_matrix[relevant_entries[i, "row"], relevant_entries[i, "col"]])
    if (organize_by == "source") {
      title <- paste0("Dest: ", dest, " (", dest_antigen, ")", "\nComp Value: ", spill_value)
    } else {
      title <- paste0("Src: ", source, " (", source_antigen, ")", "\nComp Value: ", spill_value)
    }

    p <- tryCatch({
      plotScatter(cytof_data, c(source_antigen, dest_antigen), zeros = TRUE) +
        ggtitle(title) +
        theme(plot.title   = element_text(size = title_size, hjust = 0.5),
              axis.title   = element_text(size = axis_size),
              axis.title.x = element_text(size = axis_size),
              axis.title.y = element_text(size = axis_size),
              axis.text    = element_text(size = axis_size - 1),
              axis.text.x  = element_text(size = axis_size - 1),
              axis.text.y  = element_text(size = axis_size - 1),
              plot.margin  = margin(2, 2, 2, 2, "pt"),
              legend.position = "none")
    }, error = function(e) {
      message(paste0("plotScatter failed for ", source, " to ", dest, ": ", e$message))
      NULL
    })

    if (!is.null(p) && inherits(p, "ggplot")) {
      group_key <- if (organize_by == "source") source else dest
      if (!group_key %in% names(grouped_plots)) grouped_plots[[group_key]] <- list()
      grouped_plots[[group_key]][[length(grouped_plots[[group_key]]) + 1]] <- p
    }
  }

  if (length(grouped_plots) == 0) {
    warning("No plots were generated.")
    return(NULL)
  }

  all_plots <- list()
  for (channel_name in names(grouped_plots)) {
    plist <- grouped_plots[[channel_name]]
    channel_antigen <- panel_table$antigen[panel_table$fcs_colname == channel_name]
    if (length(channel_antigen) < 1) channel_antigen <- "Unknown"
    label_text <- if (organize_by == "source") "Source:" else "Destination:"
    label_plot <- ggplot() +
      annotate("text", x = 0.5, y = 0.58, label = label_text,      size = 5, fontface = "bold", color = "darkblue") +
      annotate("text", x = 0.5, y = 0.50, label = channel_name,    size = 5, fontface = "bold", color = "darkred") +
      annotate("text", x = 0.5, y = 0.42, label = channel_antigen, size = 5, fontface = "bold", color = "darkgreen") +
      theme_void() +
      coord_cartesian(xlim = c(0, 1), ylim = c(0.35, 0.65)) +
      theme(plot.margin = margin(1, 1, 1, 1, "pt"),
            panel.background = element_rect(fill = "grey95", color = "grey80"),
            aspect.ratio = 0.8)
    all_plots <- c(all_plots, list(label_plot), plist)
  }

  cowplot::plot_grid(plotlist = all_plots, ncol = n_col, align = "hv")
}
