# plotAbundancesInhibStats
#
# Cluster-abundance strip plot with a choice of summary geom and optional
# pairwise statistics. Migrated faithfully from the nBass26 inhibitor analysis
# (250115 nBass26 inhibitor expt) — the canonical version David settled on.
#
# Key behaviour for `geom = "median_line"`: the donor points are drawn FIRST,
# then the median crossbar is layered ON TOP, so the median line sits in front
# of the points (a flat crossbar, ymin = ymax = median). `geom = "boxplot"`
# gives the older box-and-whisker style.
#
# Points are mapped by `shape_by`; add a scale_shape_manual(values = ...) on the
# returned plot to use specific (filled, 21-25) donor shapes with the grey84
# point fill.
#
# Stats (show_stats = TRUE) use rstatix (dunn_test / wilcox_test) +
# ggpubr::stat_pvalue_manual; those packages are only needed when stats are on.

plotAbundancesInhibStats <- function(x, k = "meta20", by = c("sample_id", "cluster_id"),
                               group_by = "condition", shape_by = NULL, col_clust = TRUE,
                               n_cols = 4, log = FALSE, miny = 0.01, maxy = NA,
                               point_size = 2, panel_spacing = 1, step_increase = 0.1, facet_ratio = 1, textsize = 14,
                               external_stats = NULL, show_stats = FALSE, point_alpha = 1,
                               distance = c("euclidean", "maximum", "manhattan", "canberra",
                                         "binary", "minkowski"),
                               linkage = c("average", "ward.D", "single", "complete",
                                         "mcquitty", "median", "centroid", "ward.D2"),
                               k_pal = .wl_cluster_cols, clusters_order = NULL,
                               geom = c("boxplot", "median_line"),
                               vlines = NULL) {

  library(ggh4x)
  library(rstatix)

  # Input checking
  by <- match.arg(by)
  geom <- match.arg(geom)
  .wl_check_sce(x, TRUE)
  k <- .wl_check_k(x, k)
  .wl_check_cd_factor(x, group_by)
  .wl_check_cd_factor(x, shape_by)
  .wl_check_pal(k_pal)
  linkage <- match.arg(linkage)
  distance <- match.arg(distance)
  stopifnot(is.logical(col_clust), length(col_clust) == 1)

  # Original data preparation code remains the same
  shapes <- .wl_get_shapes(x, shape_by)
  if (is.null(shapes)) shape_by <- NULL

  if (by == "sample_id") {
    nk <- nlevels(.wl_cluster_ids(x, k))
    if (length(k_pal) < nk)
      k_pal <- colorRampPalette(k_pal)(nk)
  }

  # Calculate frequencies
  ns <- table(cluster_id = .wl_cluster_ids(x, k), sample_id = .wl_sample_ids(x))
  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)
  m <- match(df$sample_id, x$sample_id)
  for (i in c(shape_by, group_by)) df[[i]] <- x[[i]][m]

  # Clustering code remains the same
  if (by == "sample_id" && col_clust && length(unique(df$sample_id)) > 1) {
    d <- dist(t(fq), distance)
    h <- hclust(d, linkage)
    o <- colnames(fq)[h$order]
    df$sample_id <- factor(df$sample_id, o)
  }

  # Apply clusters order if provided
  if (!is.null(clusters_order)) {
    df$cluster_id <- factor(df$cluster_id, levels = clusters_order)
  }

  # Log transformation
  if (log == TRUE) {
    df$Freq <- df$Freq + 0.02
  }

  dfout <<- df

  # Calculate max y values for each cluster
  maxy_values <- df %>%
    group_by(cluster_id) %>%
    dplyr::summarize(maxy = ceiling(max(Freq, na.rm = TRUE) * 1.1))

  # Create base plot
  p <- ggplot(df, aes_string(y = "Freq")) +
    labs(x = NULL, y = "Proportion [%]") +
    theme_bw() +
theme(
              panel.grid = element_blank(), # Remove grid lines
              panel.background = element_rect(fill = "transparent", color = NA), # Transparent panel background
              panel.border = element_blank(), # Remove the box around the panel
              panel.spacing = unit(panel_spacing, "lines"), # Panel spacing
              strip.background = element_rect(fill = NA, color = NA), # No background for strip
              strip.text = element_text(size = textsize, face = "bold", color = "black"), # Strip text styling
              strip.placement = "outside", # Move facet labels outside
              axis.ticks = element_line(color = "black"), # Black axis ticks
              axis.ticks.length = unit(0.3, "lines"), # Longer axis ticks
              axis.line = element_line(color = "black", size = 0.5), # Add x and y axis lines for all panels
              axis.text = element_text(color = "black", size = textsize), # Axis text styling
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), # Rotate x-axis text
              axis.text.y = element_text(color = "black"), # Y-axis text color
              axis.title = element_text(size = textsize), # Axis title styling
              legend.text = element_text(size = textsize), # Legend text styling
              legend.title = element_text(size = textsize), # Legend title styling
              plot.margin = unit(c(0, 1, 0, 0), "lines"), # Plot margin
              text = element_text(size = textsize), # General text styling
              aspect.ratio = facet_ratio # Maintain aspect ratio
)

  # Add facets
  p <- p + facet_wrap(~cluster_id, scales = "free_y", ncol = n_cols)

  # Determine fixed dodge width for consistency
  dodge_width <- 0.8

  # Add points first (hacked because we are just using median_line for now)
  if (!is.null(shape_by)) {
    if (geom == "median_line") {
      p <- p + geom_quasirandom(
        aes_string(x = group_by, shape = shape_by, group = group_by),
        fill = "grey84",
        size = point_size,
        width = 0.2,
        dodge.width = dodge_width,
        alpha = point_alpha
      )
    } else {
      p <- p + geom_quasirandom(
        aes_string(x = group_by, shape = shape_by),
        fill = "grey84",
        size = point_size,
        width = 0.2,
        alpha = point_alpha
      )
    }
  } else {
    if (geom == "median_line") {
      p <- p + geom_quasirandom(
        aes_string(x = group_by, group = group_by),
        fill = "grey84",
        size = point_size,
        width = 0.2,
        shape = 21,
        dodge.width = dodge_width,
        alpha = point_alpha
      )
    } else {
      p <- p + geom_quasirandom(
        aes_string(x = group_by),
        fill = "grey84",
        size = point_size,
        width = 0.2,
        shape = 21,
        alpha = point_alpha
      )
    }
  }

  # Add boxplot or median line based on geom parameter
  if (geom == "boxplot") {
    p <- p + geom_boxplot(aes_string(x = group_by, fill = group_by),
                         color = "black", position = position_dodge(),
                         alpha = 0.8, outlier.color = NA, show.legend = FALSE)
  } else if (geom == "median_line") {
    # Create a grouping variable for median lines
    df$group_median <- interaction(df[[group_by]], df$cluster_id)

    # Add median crossbar
    p <- p + stat_summary(
      aes_string(x = group_by, group = group_by, fill = group_by),
      fun.data = function(x) {
        m <- median(x, na.rm = TRUE)
        data.frame(y = m, ymin = m, ymax = m)
      },
      geom = "crossbar",
      width = 0.7,
      color = "black",
      size = 0.8,
      fatten = 1,
      position = position_dodge(width = dodge_width),
      alpha = 0.8,
      # [ggplot2 4.0] the crossbar/boxplot key otherwise renders INTO the shape
      # (e.g. patient_id) legend and occludes the point glyphs. These are visual
      # indicators that don't need their own legend entry.
      show.legend = FALSE
    )
  }

  # Add vertical lines between conditions at specified positions
  if (!is.null(vlines)) {
    # Get unique condition values
    conditions <- sort(unique(df[[group_by]]))

    # Create a data frame for vertical lines
    vline_positions <- data.frame()

    # Process each vline specification (position between which conditions to draw lines)
    for (pos in vlines) {
      if (pos > 0 && pos < length(conditions)) {
        # Calculate position between the two conditions (pos and pos+1)
        x_val <- pos + 0.5
        vline_positions <- rbind(vline_positions, data.frame(x = x_val))
      }
    }

    if (nrow(vline_positions) > 0) {
      # Add vertical lines that span from top to bottom of each facet
      p <- p + geom_vline(
        data = vline_positions,
        aes(xintercept = x),
        linetype = "dashed",
        color = "darkgray",
        size = 0.8
      )
    }
  }

  # Add statistics if requested
  if (show_stats) {
    if (!is.null(external_stats)) {
      # Use external stats if provided
      dummy_stat <- df %>%
        group_by(cluster_id) %>%
        wilcox_test(as.formula(paste("Freq ~", group_by)), paired = FALSE)
      dummy_stat <- dummy_stat %>%
        add_y_position(scales = "free", step.increase = step_increase)

      # Merge y.position from dummy_stat to external_stats
      external_stats <- external_stats %>%
        left_join(dummy_stat %>% select(cluster_id, group1, group2, y.position),
                 by = c("cluster_id", "group1", "group2"))

      stat.test <- external_stats
    } else {
      # Calculate stats internally
      stat.test <- df %>%
        group_by(cluster_id) %>%
        dunn_test(as.formula(paste("Freq ~", group_by)), p.adjust.method = "holm") %>%
        add_y_position(scales = "free", step.increase = step_increase)
    }

    statout <<- stat.test
    p <- p + stat_pvalue_manual(stat.test, label = "p = {scales::pvalue(p.adj)}",
                               tip.length = 0.01, hide.ns = FALSE, size = 4)
  }

  # Add scale transformations
  if (log == TRUE) {
    p <- p +
      scale_y_continuous(trans = 'log10',
                        limits = c(miny, maxy),
                        breaks = c(0.01, 0.1, 1, 10, 100),
                        labels = c(0.01, 0.1, 1, 10, 100)) +
      annotation_logticks(base = 10, sides = "l", outside = TRUE) +
      coord_cartesian(clip = "off") +
      scale_size_area(max_size = 15) +
      theme(axis.text.y = element_text(margin = margin(r = 8)))
  } else {
    p <- p +
      scale_y_continuous(limits = c(0, maxy)) +
      coord_cartesian(clip = "off") +
      scale_size_area(max_size = 15)
  }

  return(p)
}

