# plotIsotypeByConditionDay — migrated from CyTOF nBass_helpers.R into seekit (CATALYST-free).
# 2026-06-10: lifted verbatim, de-CATALYST'd (.wl_* internals, namespace hack removed).

plotIsotypeByConditionDay <- function(sce,
                                     master_cluster_column,
                                     master_cluster_value,
                                     daughter_cluster_column = "isotype",
                                     condition_day_column = "condition_day",
                                     sample_id_column = "sample_id",
                                     patient_id_column = "patient_id",
                                     textsize = 18,
                                     panel_spacing = 1,
                                     facet_ratio = 1,
                                     fill_by = "patient_id",
                                     shape_by = NULL,
                                     point_size = 3,
                                     point_alpha = 0.9,
                                     cell_threshold = 100,
                                     condition_day_order = NULL,
                                     daughter_order = NULL,
                                     n_cols = 4) {

  df <- as.data.frame(colData(sce))

  if (!is.null(master_cluster_value)) {
    df <- df[df[[master_cluster_column]] == master_cluster_value, ]
    if(nrow(df) == 0) {
      return(ggplot() + theme_void() + ggtitle(paste("No data available for", master_cluster_value)))
    }
  }

  if (!is.null(daughter_order)) {
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = daughter_order)
  } else {
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]])
  }

  if (!is.null(condition_day_order)) {
    if (is.data.frame(condition_day_order) || is.list(condition_day_order)) {
      condition_levels <- as.character(unlist(condition_day_order))
    } else {
      condition_levels <- as.character(condition_day_order)
    }
    df[[condition_day_column]] <- factor(df[[condition_day_column]], levels = condition_levels)
  } else {
    df[[condition_day_column]] <- factor(df[[condition_day_column]])
  }

  cluster_ids_abundances <- df[[daughter_cluster_column]]

  ns <- table(cluster_id = cluster_ids_abundances, sample_id = df[[sample_id_column]])

  df_abund <- as.data.frame(ns)
  names(df_abund) <- c(daughter_cluster_column, sample_id_column, "cell_count")

  meta_info <- unique(df[, c(sample_id_column, patient_id_column, condition_day_column)])
  df_abund <- merge(df_abund, meta_info, by = sample_id_column)

  sample_totals <- df_abund %>%
    dplyr::group_by(dplyr::across(c(dplyr::all_of(c(sample_id_column, condition_day_column))))) %>%
    dplyr::summarise(total_cells = sum(cell_count), .groups = "drop") %>%
    dplyr::mutate(keep = total_cells >= cell_threshold)

  excluded <- sample_totals %>% dplyr::filter(!keep)
  if(nrow(excluded) > 0) {
    message("Excluded samples (threshold=", cell_threshold, "):")
    for(i in 1:nrow(excluded)) {
      message("Sample: ", excluded[[sample_id_column]][i],
              ", ", condition_day_column, ": ", excluded[[condition_day_column]][i],
              ", Total cells: ", excluded$total_cells[i])
    }
  }

  valid_samples <- sample_totals %>%
    dplyr::filter(keep) %>%
    dplyr::select(dplyr::all_of(c(sample_id_column, condition_day_column)))

  df_filtered <- df_abund %>%
    dplyr::inner_join(valid_samples, by = c(sample_id_column, condition_day_column))

  if(nrow(df_filtered) == 0) {
    return(ggplot() + theme_void() +
             ggtitle(paste("No samples meet the threshold of", cell_threshold, "cells")))
  }

  plot_df <- df_filtered %>%
    dplyr::group_by(dplyr::across(c(dplyr::all_of(c(sample_id_column, condition_day_column))))) %>%
    dplyr::mutate(proportion = cell_count / sum(cell_count) * 100) %>%
    dplyr::ungroup()

  median_df <- plot_df %>%
    dplyr::group_by(dplyr::across(c(dplyr::all_of(c(daughter_cluster_column, condition_day_column))))) %>%
    dplyr::summarise(
      median_prop = median(proportion, na.rm = TRUE),
      .groups = "drop"
    )

  median_df[[condition_day_column]] <- factor(median_df[[condition_day_column]],
                                              levels = levels(df[[condition_day_column]]))

  plot_df[[condition_day_column]] <- factor(plot_df[[condition_day_column]],
                                            levels = levels(df[[condition_day_column]]))
  plot_df[[daughter_cluster_column]] <- factor(plot_df[[daughter_cluster_column]],
                                               levels = levels(df[[daughter_cluster_column]]))

  dodge_width <- 0.5

  p <- ggplot() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      panel.border = element_blank(),
      panel.spacing = unit(panel_spacing, "lines"),
      strip.background = element_rect(fill = NA, color = NA),
      strip.text = element_text(size = textsize, face = "bold", color = "black"),
      strip.placement = "outside",
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(0.3, "lines"),
      axis.line = element_line(color = "black", size = 0.5),
      axis.text = element_text(color = "black", size = textsize),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(color = "black"),
      axis.title = element_text(size = textsize),
      legend.text = element_text(size = textsize),
      legend.title = element_text(size = textsize),
      plot.margin = unit(c(0, 1, 0, 0), "lines"),
      text = element_text(size = textsize),
      aspect.ratio = facet_ratio
    )

  if (!is.null(shape_by)) {
    p <- p + ggbeeswarm::geom_quasirandom(
      data = plot_df,
      aes(x = .data[[condition_day_column]],
          y = proportion,
          shape = .data[[shape_by]],
          group = interaction(.data[[condition_day_column]], .data[[fill_by]])),
      color = "black",
      fill = "grey84",
      stroke = 0.5,
      dodge.width = dodge_width,
      width = 0.1,
      groupOnX = TRUE,
      size = point_size,
      alpha = point_alpha,
      show.legend = TRUE
    ) +
      scale_shape_manual(values = c(21, 22, 24, 25))
  } else {
    p <- p + ggbeeswarm::geom_quasirandom(
      data = plot_df,
      aes(x = .data[[condition_day_column]],
          y = proportion,
          fill = .data[[fill_by]],
          group = .data[[condition_day_column]]),
      shape = 21,
      color = "black",
      stroke = 0.5,
      width = 0.1,
      size = point_size,
      alpha = point_alpha
    )
  }

  p <- p + geom_crossbar(
    data = median_df,
    aes(x = .data[[condition_day_column]],
        y = median_prop,
        ymin = median_prop,
        ymax = median_prop),
    width = 0.7,
    color = "black",
    size = 0.8,
    fatten = 1
  )

  p <- p +
    # axes = "y": keep a y-axis on every panel (needed for free_y), but draw the
    # x-axis (long condition_day labels) only on the bottom row so it isn't squashed.
    facet_wrap2(~ .data[[daughter_cluster_column]], scales = "free_y", ncol = n_cols, axes = "y") +
    labs(
      x = condition_day_column,
      y = "Proportion [%]",
      fill = fill_by,
      title = ifelse(!is.null(master_cluster_value),
                     paste("Isotype distribution in", master_cluster_column, "=", master_cluster_value),
                     "Isotype distribution across conditions")
    ) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)))

  p$data <- plot_df
  return(p)
}
