# plotDR2.R — CATALYST-free rewrite of plotDR2, with plotDR1's unique features
# folded in (border_width, plot_order, highlight_cluster) so plotDR1 can retire to
# a deprecated alias -> plotDR2.
# -----------------------------------------------------------------------------
# Body is David Priest's own code (merged from plotDR1 + plotDR2). Only the
# CATALYST internal calls are swapped for the MIT .wl_* equivalents; ggplot2/
# reshape2/ggrastr stay bare (imported by the package). No CATALYST dependency.
#
# Superset arg map:
#   from plotDR2 : parse_color_by, legend_title, cluster_order, colData clusterings
#   from plotDR1 : border_width, plot_order, highlight_cluster
# -----------------------------------------------------------------------------

plotDR2 <- function(x, dr = NULL, textsize = 18, legendpointsize = 7,
                       color_by = "condition", facet_by = NULL, hide_axis = FALSE,
                       alpha = 0.8, pointsize = 0.4, ncol = NULL, assay = "exprs",
                       scale = TRUE, random_order = FALSE, q = 0.01, dims = c(1, 2),
                       border_width = 1, k_pal = .wl_cluster_cols,
                       a_pal = hcl.colors(10, "Viridis"), rast = FALSE,
                       panel_spacing = 1, parse_color_by = FALSE, legend_title = NULL,
                       cluster_order = NULL, plot_order = NULL, highlight_cluster = NULL) {

  stopifnot(is(x, "SingleCellExperiment"), .wl_check_assay(x, assay),
            length(reducedDims(x)) != 0, is.logical(scale), length(scale) == 1,
            is.numeric(q), length(q) == 1, q >= 0, q < 0.5)

  .wl_check_pal(a_pal)
  .wl_check_cd_factor(x, facet_by)

  if (!is.null(ncol)) stopifnot(is.numeric(ncol), length(ncol) == 1, ncol %% 1 == 0)

  if (is.null(dr)) {
    dr <- reducedDimNames(x)[1]
  } else {
    stopifnot(is.character(dr), length(dr) == 1, dr %in% reducedDimNames(x))
  }

  stopifnot(is.numeric(dims), length(dims) == 2, dims %in% seq_len(ncol(reducedDim(x, dr))))

  if (!all(color_by %in% rownames(x))) {
    stopifnot(length(color_by) == 1)
    if (!color_by %in% names(colData(x))) {
      .wl_check_sce(x, TRUE)
      .wl_check_pal(k_pal)
      .wl_check_k(x, color_by)
      kids <- .wl_cluster_ids(x, color_by)
      if (!is.null(cluster_order)) kids <- factor(kids, levels = cluster_order)
      nk <- nlevels(kids)
      if (length(k_pal) < nk) k_pal <- colorRampPalette(k_pal)(nk)
      # [#3] dropped the `plotDR1_colorpal <<-` global side-effect — use
      #      cluster_palette()/cluster_palette() to get a named palette instead.
    } else {
      kids <- factor(colData(x)[[color_by]])
      if (!is.null(cluster_order)) kids <- factor(kids, levels = cluster_order)
    }
  }

  xy <- reducedDim(x, dr)[, dims]
  colnames(xy) <- c("x", "y")
  df <- data.frame(colData(x), xy, check.names = FALSE)

  if (all(color_by %in% rownames(x))) {
    es <- as.matrix(assay(x, assay))
    es <- es[color_by, , drop = FALSE]
    if (scale) es <- .wl_scale_exprs(es, 1, q)
    df <- reshape2::melt(cbind(df, t(es)), id.vars = colnames(df))
    l <- switch(assay, exprs = "expression", assay)
    l <- paste0("scaled\n"[scale], l)
    scale <- scale_colour_gradientn(l, colors = a_pal)
    thm <- guide <- NULL
    color_by <- "value"
    facet <- facet_wrap("variable", ncol = ncol)
  } else if (is.numeric(df[[color_by]])) {
    if (scale) {
      vs <- as.matrix(df[[color_by]])
      df[[color_by]] <- .wl_scale_exprs(vs, 2, q)
    }
    l <- paste0("scaled\n"[scale], color_by)
    scale <- scale_colour_gradientn(l, colors = a_pal)
    # [#4] no backtick-quoting needed: aes(.data[[color_by]]) handles any name
    facet <- thm <- guide <- NULL
  } else {
    facet <- NULL
    if (exists("kids")) {
      df[[color_by]] <- kids
      scale <- scale_color_manual(values = k_pal)
    } else scale <- NULL

    n <- nlevels(droplevels(factor(df[[color_by]])))
    guide <- guides(col = guide_legend(ncol = ifelse(n > 12, 2, 1),
                                       override.aes = list(alpha = 1, size = legendpointsize)))
    thm <- theme(legend.key.height = unit(0.8, "lines"), text = element_text(size = textsize))
  }

  if (dr %in% c("PCA", "MDS")) {
    asp <- coord_equal()
  } else asp <- NULL

  if (dr == "PCA") {
    labs <- paste0("PC", dims)
  } else labs <- paste(dr, "dim.", dims)

  df <- df[!(is.na(df$x) | is.na(df$y)), ]

  # Point draw order: explicit factor order (plot_order) takes precedence,
  # else optional shuffle (random_order). [from plotDR1]
  if (!is.null(plot_order)) {
    df[[color_by]] <- factor(df[[color_by]], levels = plot_order)
    df <- df[order(df[[color_by]]), ]
  } else if (random_order) {
    set.seed(1)
    df <- df[sample(nrow(df)), ]
  }

  # [#4] aes_string (deprecated) -> tidy-eval aes(.data[[...]])
  p <- ggplot(df, aes(x = .data[["x"]], y = .data[["y"]], colour = .data[[color_by]]))

  if (rast) {
    # Cairo-FREE rasterisation via scattermore (.wl_geom_point_rast) — never uses
    # ggrastr/Cairo, which fails on large UMAPs.
    p <- p + .wl_geom_point_rast(size = pointsize, alpha = alpha, shape = 16, raster.dpi = 600)
  } else {
    p <- p + geom_point(size = pointsize, alpha = alpha, shape = 16)
  }

  # Highlight one cluster with outlined points + density contours [from plotDR1]
  if (!is.null(highlight_cluster)) {
    p <- p + geom_point(
        data = df[df[[color_by]] == highlight_cluster, ],
        aes(x = x, y = y),
        size = 2, alpha = 0.8, shape = 21, color = "black", fill = "yellow") +
      stat_density_2d(
        data = df[df[[color_by]] == highlight_cluster, ],
        aes(x = x, y = y, fill = after_stat(level)),
        geom = "polygon", color = "black", alpha = 0.3)
  }

  p <- p + labs(x = labs[1], y = labs[2]) + facet + scale + guide + asp +
    theme_minimal() + thm +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.grid.major = element_blank(),
          axis.text = element_text(color = "black"),
          panel.spacing = unit(panel_spacing, "lines"),
          aspect.ratio = if (is.null(asp)) 1 else NULL)

  # Panel border width is now tunable [from plotDR1]
  p <- p + theme(panel.border = element_rect(colour = "black", fill = NA,
                                              linewidth = border_width))

  if (hide_axis == TRUE) {
    p <- p + theme(axis.text.x = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.ticks.y = element_blank())
  }

  if (!is.null(legend_title)) {
    p <- p + guides(col = guide_legend(title = legend_title,
                                       override.aes = list(alpha = 1, size = legendpointsize)))
  }

  if (parse_color_by) {
    p <- p + scale_color_manual(values = k_pal,
      labels = sapply(levels(df[[color_by]]), function(name) as.expression(parse(text = name))))
  }

  if (is.null(facet_by)) return(p)

  if (is.null(facet)) {
    p + facet_wrap(facet_by, ncol = ncol)
  } else {
    if (nlevels(df$variable) == 1) {
      p + facet_wrap(facet_by, ncol = ncol) + ggtitle(levels(df$variable))
    } else {
      fs <- c("variable", facet_by)
      ns <- vapply(df[fs], nlevels, numeric(1))
      if (ns[2] > ns[1]) fs <- rev(fs)
      p + facet_grid(reformulate(fs[1], fs[2]))
    }
  }
}
