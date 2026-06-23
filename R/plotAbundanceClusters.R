# plotAbundanceClusters.R — CATALYST-free.
#
# Per-cluster abundance over a time/grouping variable (e.g. day), one facet per
# cluster, points coloured by a condition variable with an optional connecting
# line and mean +/- SEM error bars. Migrated from the iGCB nBass CyTOF analysis
# (241119 nBass Analysis 2); used for the cluster-proportion-over-days panels
# (Fig 1 proportions / supplements).
#
# The function BODY is David Priest's own code; it is not derived from CATALYST.
# Its calls into CATALYST internals (.check_sce/.check_k/.check_cd_factor) and
# the CATALYST generics cluster_ids()/sample_ids() are replaced with the MIT
# .wl_* equivalents (R/wl_internals.R), and the
#   environment(...) <- asNamespace('CATALYST')
# hack is removed. All other symbols resolve via the package NAMESPACE imports
# (ggplot2, dplyr, rlang, ggh4x, ggbeeswarm, grid), so the function needs no
# CATALYST on the search path.
#
# geom = "boxes" | "bar" | "none" sets the per-day summary geom; geom_type =
# "quasirandom" | "point" sets the per-sample points (drawn on top). average_samples
# = TRUE collapses to per-(condition,day,cluster) means and adds SEM error bars.

plotAbundanceClusters <- function(x,
                                  k = "somPB",
                                  selected_cluster = "PB",  # This parameter is now ignored.
                                  day_var = "day",
                                  color_var = "condition",
                                  point_size = 2,
                                  point_stroke = 1,
                                  panel_spacing = 1,
                                  textsize = 16,
                                  clusters_order = NULL,
                                  log = FALSE,
                                  n_cols = 4,
                                  miny = 0.01,
                                  maxy = NA,
                                  average_samples = FALSE,
                                  geom = c("boxes", "bar", "none"),
                                  mean_or_med = "mean", # This only says where box or barplot should be (not about averaged points)
                                  facet_ratio = 1.5,
                                  line = FALSE,
                                  geom_type = c("quasirandom", "point")) {
  # Check input
  .wl_check_sce(x, TRUE)
  k <- .wl_check_k(x, k)
  .wl_check_cd_factor(x, day_var)
  .wl_check_cd_factor(x, color_var)

  # Compute the cluster frequencies per sample
  ns <- table(cluster_id = .wl_cluster_ids(x, k), sample_id = .wl_sample_ids(x))
  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)

  # Add day and condition metadata for each sample
  m <- match(df$sample_id, x$sample_id)
  df[[day_var]] <- x[[day_var]][m]
  df[[color_var]] <- x[[color_var]][m]

  if (log) {
    df$Freq <- df$Freq + 0.02
  }

  # For storing original data (needed for error bars)
  df_original <- df

  # Average across samples if requested
  if (average_samples) {
    # Calculate summary statistics for error bars
    df_summary <- df_original %>%
      dplyr::group_by(across(all_of(c(color_var, day_var, "cluster_id")))) %>%
      dplyr::summarise(
        mean_freq = mean(Freq, na.rm = TRUE),
        sem = stats::sd(Freq, na.rm = TRUE) / sqrt(sum(!is.na(Freq))),
        .groups = "drop"
      )

    # Replace main dataframe with averaged values
    df <- df_original %>%
      dplyr::group_by(across(all_of(c(color_var, day_var, "cluster_id")))) %>%
      dplyr::summarise(Freq = mean(Freq, na.rm = TRUE), .groups = "drop")
  }

  geom <- match.arg(geom)
  mean_or_med <- match.arg(mean_or_med)
  fun <- if (mean_or_med == "median") "median" else "mean"
  geom_type <- match.arg(geom_type)

  if (!is.null(clusters_order)) {
    df$cluster_id <- factor(df$cluster_id, levels = clusters_order, ordered = TRUE)
    if (average_samples) {
      df_summary$cluster_id <- factor(df_summary$cluster_id, levels = clusters_order, ordered = TRUE)
    }
  }

  # Create the plot with facet_wrap2 (ggh4x) to show x-axis ticks on all panels.
  p <- ggplot(df, aes_string(x = day_var, y = "Freq", color = color_var)) +
    labs(x = day_var, y = "Proportion [%]") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          aspect.ratio = facet_ratio,
          axis.text = element_text(color = "black", size = textsize),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
    ggh4x::facet_wrap2(~ cluster_id, scales = "free_y", axes = "all", ncol = n_cols)

  # Build plot incrementally to control the order of layers

  # First add boxes or bars based on geom selection - now supporting "none"
  if (geom == "boxes") {
    p <- p + geom_boxplot(aes_string(group = day_var),
                          fill = "grey84", color = "black",
                          width = 0.75, linewidth = 0.3, alpha = 0.9,
                          outlier.color = NA, show.legend = FALSE)
  } else if (geom == "bar") {
    p <- p + stat_summary(fun = fun, geom = "bar",
                          aes_string(group = day_var),
                          fill = "grey84", color = "black",
                          position = position_dodge(), alpha = 1)
  }
  # 'none' option doesn't add any geom

  # Add line if requested with grey colour
  if (line) {
    p <- p + geom_line(aes_string(group = color_var), color = "grey",
                       linewidth = 0.4, alpha = 0.6)
  }

  # Add points based on geom_type selection (now added last so they're on top)
  if (geom_type == "quasirandom") {
    p <- p + ggbeeswarm::geom_quasirandom(aes_string(fill = color_var, color = color_var),
                              size = point_size, width = 0.4,
                              stroke = point_stroke, show.legend = TRUE, shape = 21)
  } else if (geom_type == "point") {
    p <- p + geom_point(aes_string(fill = color_var, color = color_var),
                         size = point_size, stroke = point_stroke,
                        show.legend = TRUE, shape = 21)
  }

  # Add error bars when averaging samples (placed before points so they appear behind)
  if (average_samples) {
    # Add error bars as a separate layer - now with black color
    p <- p + geom_errorbar(
      data = df_summary,
      aes(x = .data[[day_var]],
          y = mean_freq,
          ymin = mean_freq - sem,
          ymax = mean_freq + sem,
          group = .data[[color_var]]),
      color = "black", # Fixed black color for all error bars
      width = 0.2,
      linewidth = 0.5
    )
  }

  # Apply log-scale or linear scale with expansion around y=0 and top
  if (log) {
    p <- p + scale_y_continuous(trans = 'log10', limits = c(miny, maxy),
                                breaks = c(0.01, 0.1, 1, 10, 100),
                                labels = c(0.01, 0.1, 1, 10, 100),
                                expand = expansion(mult = c(0.1, 0.1))) +
      annotation_logticks(sides = "l", outside = TRUE) +
      coord_cartesian(clip = "off")
  } else {
    p <- p + scale_y_continuous(limits = c(0, maxy),
                                expand = expansion(mult = c(0.1, 0.1))) +
      coord_cartesian(clip = "off")
  }

p <- p + theme(
              strip.background = element_blank(),
              strip.text = element_text(face = "bold", color = "black", size = textsize),
              panel.spacing = unit(panel_spacing, "lines"),
              panel.background = element_rect(fill = "transparent", color = NA),
              axis.ticks = element_line(color = "black"), # Ensure ticks are black for both axes
              axis.ticks.length = unit(0.3, "lines"), # Make the ticks longer
              axis.text.x = element_text(color = "black"),
              plot.margin = unit(c(0, 1, 0, 0), "lines"),
              axis.line.x = element_line(color = "black", size = 0.4), # Add x-axis line
              axis.line.y = element_line(color = "black", size = 0.4), # Add y-axis line
              panel.border = element_blank(), # Remove the box around the panel
              axis.text.y = element_text(color = "black"),
              text = element_text(size = textsize))

return(p)

}
