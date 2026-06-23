#' @title FeaturePlot2
#' @description Seurat::FeaturePlot drop-in with two David-flavoured
#'   additions on top of stock Seurat 5.x:
#'   - `myTheme = TRUE` (default): when there is no `split.by`, applies the
#'     standard publication theme used across the lab -- 1:1 aspect ratio,
#'     panel border, no axis ticks/labels/lines, adjustable `textsize`.
#'   - `alpha` is honoured on the points.
#'
#'   All other arguments mirror Seurat::FeaturePlot. Function environment
#'   is rebound to `asNamespace("Seurat")` at the bottom of the file so
#'   the internals it calls (SingleDimPlot, CenterTitle, LabelClusters,
#'   SetQuantile, BlendMatrix, BlendMap, BlendExpression, RandomName,
#'   NoLegend, IFeaturePlot, `%||%`, `%iff%`, `brewer.pal.info`, ...)
#'   resolve without explicit `:::`.
#'
#' @param object Seurat object.
#' @param features Features to plot (genes, metadata, dim-reduction columns).
#' @param dims Two-length integer of reduction dims. Default: c(1, 2)
#' @param cells Cells to include. Default: NULL
#' @param cols Colour scale. Default: lightgrey -> blue (lightgrey/red/green if blend).
#' @param pt.size Point size. Default: NULL
#' @param alpha Point alpha. Default: 1
#' @param order Plot high-expression cells last. Default: FALSE
#' @param min.cutoff Per-feature minimum cutoff (NA = observed min). Default: NA
#' @param textsize Base text size used by myTheme branch. Default: 12
#' @param max.cutoff Per-feature maximum cutoff (NA = observed max). Default: NA
#' @param reduction Reduction name. Default: NULL (Seurat default)
#' @param split.by Metadata column to split by. Default: NULL
#' @param keep.scale "feature" or "all"; controls scale sharing. Default: "feature"
#' @param shape.by Metadata column for point shape. Default: NULL
#' @param slot Assay slot to pull from. Default: "data"
#' @param blend Two-feature blend plot. Default: FALSE
#' @param blend.threshold Blend threshold. Default: 0.5
#' @param myTheme Apply the lab publication theme (border, no axes,
#'   aspect.ratio = 1). Only takes effect when there is no split.by. Default: TRUE
#' @param label Add cluster labels. Default: FALSE
#' @param label.size Label text size. Default: 4
#' @param label.color Label colour. Default: "black"
#' @param repel Use ggrepel for labels. Default: FALSE
#' @param ncol Patchwork columns. Default: NULL (auto)
#' @param coord.fixed Force coord_fixed(). Default: FALSE (myTheme already locks aspect)
#' @param by.col Split layout direction. Default: TRUE
#' @param sort.cell Deprecated -- use `order`.
#' @param interactive Open IFeaturePlot. Default: FALSE
#' @param combine Combine via patchwork. Default: TRUE
#' @param raster Rasterise points. Default: NULL
#' @param raster.dpi Raster dpi as c(width, height). Default: c(512, 512)
#'
#' @return A ggplot (or patchwork) object.
#' @export
FeaturePlot2 <- function(object, features, dims = c(1, 2), cells = NULL,
                         cols = if (blend) {
                           c("lightgrey", "#ff0000", "#00ff00")
                         } else {
                           c("lightgrey", "blue")
                         },
                         pt.size = NULL, alpha = 1, order = FALSE,
                         min.cutoff = NA, textsize = 12, max.cutoff = NA,
                         reduction = NULL, split.by = NULL,
                         keep.scale = "feature", shape.by = NULL,
                         slot = "data", blend = FALSE, blend.threshold = 0.5,
                         myTheme = TRUE, label = FALSE, label.size = 4,
                         label.color = "black", repel = FALSE, ncol = NULL,
                         coord.fixed = FALSE, by.col = TRUE,
                         sort.cell = deprecated(), interactive = FALSE,
                         combine = TRUE, raster = NULL,
                         raster.dpi = c(512, 512)) {
  if (is_present(arg = sort.cell)) {
    deprecate_stop(when = "4.9.0", what = "FeaturePlot(sort.cell = )",
                   with = "FeaturePlot(order = )")
  }
  if (isTRUE(x = interactive)) {
    return(IFeaturePlot(object = object, feature = features[1],
                        dims = dims, reduction = reduction, slot = slot))
  }
  if (!is.null(x = keep.scale)) {
    keep.scale <- arg_match0(arg = keep.scale, values = c("feature", "all"))
  }
  no.right <- theme(axis.line.y.right = element_blank(),
                    axis.ticks.y.right = element_blank(),
                    axis.text.y.right = element_blank(),
                    axis.title.y.right = element_text(face = "bold", size = 14,
                                                      margin = margin(r = 7)))
  reduction <- reduction %||% DefaultDimReduc(object = object)
  if (!is_integerish(x = dims, n = 2L, finite = TRUE) && !all(dims > 0L)) {
    abort(message = "'dims' must be a two-length integer vector")
  }
  if (isTRUE(x = blend) && length(x = features) != 2) {
    abort(message = "Blending feature plots only works with two features")
  }
  if (isTRUE(x = blend)) {
    default.colors <- eval(expr = formals(fun = FeaturePlot)$cols)
    cols <- switch(EXPR = as.character(x = length(x = cols)),
                   `0` = { warn(message = "No colors provided, using default colors"); default.colors },
                   `1` = { warn(message = paste("Only one color provided, assuming",
                                                sQuote(x = cols),
                                                "is double-negative and augmenting with default colors"))
                           c(cols, default.colors[2:3]) },
                   `2` = { warn(message = paste("Only two colors provided, assuming specified are for features and agumenting with",
                                                sQuote(default.colors[1]), "for double-negatives"))
                           c(default.colors[1], cols) },
                   `3` = cols,
                   { warn(message = "More than three colors provided, using only first three")
                     cols[1:3] })
  }
  if (isTRUE(x = blend) && length(x = cols) != 3) {
    abort("Blending feature plots only works with three colors; first one for negative cells")
  }
  dims <- paste0(Key(object = object[[reduction]]), dims)
  cells <- cells %||% Cells(x = object[[reduction]])
  data <- FetchData(object = object, vars = c(dims, "ident", features),
                    cells = cells, layer = slot)
  if (ncol(x = data) < 4) {
    abort(message = paste("None of the requested features were found:",
                          paste(features, collapse = ", "), "in slot ", slot))
  } else if (!all(dims %in% colnames(x = data))) {
    abort(message = "The dimensions requested were not found")
  }
  features <- setdiff(x = names(x = data), y = c(dims, "ident"))
  min.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = min(data[, feature]), no = cutoff))
  }, cutoff = min.cutoff, feature = features)
  max.cutoff <- mapply(FUN = function(cutoff, feature) {
    return(ifelse(test = is.na(x = cutoff), yes = max(data[, feature]), no = cutoff))
  }, cutoff = max.cutoff, feature = features)
  check.lengths <- unique(x = vapply(X = list(features, min.cutoff, max.cutoff),
                                     FUN = length, FUN.VALUE = numeric(length = 1)))
  if (length(x = check.lengths) != 1) {
    abort(message = "There must be the same number of minimum and maximum cuttoffs as there are features")
  }
  names(x = min.cutoff) <- names(x = max.cutoff) <- features
  brewer.gran <- ifelse(test = length(x = cols) == 1,
                        yes = brewer.pal.info[cols, ]$maxcolors,
                        no = length(x = cols))
  for (i in seq_along(along.with = features)) {
    f <- features[i]
    data.feature <- data[[f]]
    min.use <- SetQuantile(cutoff = min.cutoff[f], data = data.feature)
    max.use <- SetQuantile(cutoff = max.cutoff[f], data = data.feature)
    data.feature[data.feature < min.use] <- min.use
    data.feature[data.feature > max.use] <- max.use
    if (brewer.gran != 2) {
      data.feature <- if (all(data.feature == 0)) {
        rep_len(x = 0, length.out = length(x = data.feature))
      } else {
        as.numeric(x = as.factor(x = cut(x = as.numeric(x = data.feature), breaks = 2)))
      }
    }
    data[[f]] <- data.feature
  }
  data$split <- if (is.null(x = split.by)) {
    RandomName()
  } else {
    switch(EXPR = split.by,
           ident = Idents(object = object)[cells, drop = TRUE],
           object[[split.by, drop = TRUE]][cells, drop = TRUE])
  }
  if (!is.factor(x = data$split)) {
    data$split <- factor(x = data$split)
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  plots <- vector(mode = "list",
                  length = ifelse(test = blend, yes = 4,
                                  no = length(x = features) * length(x = levels(x = data$split))))
  xlims <- c(floor(x = min(data[, dims[1]])), ceiling(x = max(data[, dims[1]])))
  ylims <- c(floor(min(data[, dims[2]])), ceiling(x = max(data[, dims[2]])))
  if (blend) {
    ncol <- 4
    color.matrix <- BlendMatrix(two.colors = cols[2:3],
                                col.threshold = blend.threshold,
                                negative.color = cols[1])
    cols <- cols[2:3]
    colors <- list(color.matrix[, 1], color.matrix[1, ], as.vector(x = color.matrix))
  }
  for (i in 1:length(x = levels(x = data$split))) {
    ident <- levels(x = data$split)[i]
    data.plot <- data[as.character(x = data$split) == ident, , drop = FALSE]
    if (isTRUE(x = blend)) {
      features <- features[1:2]
      no.expression <- features[colMeans(x = data.plot[, features]) == 0]
      if (length(x = no.expression) != 0) {
        abort(message = paste("The following features have no value:",
                              paste(no.expression, collapse = ", ")))
      }
      data.plot <- cbind(data.plot[, c(dims, "ident")],
                         BlendExpression(data = data.plot[, features[1:2]]))
      features <- colnames(x = data.plot)[4:ncol(x = data.plot)]
    }
    for (j in 1:length(x = features)) {
      feature <- features[j]
      if (isTRUE(x = blend)) {
        cols.use <- as.numeric(x = as.character(x = data.plot[, feature])) + 1
        cols.use <- colors[[j]][sort(x = unique(x = cols.use))]
      } else {
        cols.use <- NULL
      }
      data.single <- data.plot[, c(dims, "ident", feature, shape.by)]
      plot <- SingleDimPlot(data = data.single, dims = dims, col.by = feature,
                            order = order, pt.size = pt.size, alpha = alpha,
                            cols = cols.use, shape.by = shape.by, label = FALSE,
                            raster = raster, raster.dpi = raster.dpi) +
        scale_x_continuous(limits = xlims) +
        scale_y_continuous(limits = ylims) +
        theme_cowplot() + CenterTitle()
      if (isTRUE(x = label)) {
        plot <- LabelClusters(plot = plot, id = "ident", repel = repel,
                              size = label.size, color = label.color)
      }
      if (length(x = levels(x = data$split)) > 1) {
        plot <- plot + theme(panel.border = element_rect(fill = NA, colour = "black"),
                             linewidth = 1)
        plot <- plot + if (i == 1) { labs(title = feature) } else { labs(title = NULL) }
        if (j == length(x = features) && !blend) {
          suppressMessages(expr = plot <- plot +
                             scale_y_continuous(sec.axis = dup_axis(name = ident),
                                                limits = ylims) + no.right)
        }
        if (j != 1) {
          plot <- plot + theme(axis.line.y = element_blank(),
                               axis.ticks.y = element_blank(),
                               axis.text.y = element_blank(),
                               axis.title.y.left = element_blank())
        }
        if (i != length(x = levels(x = data$split))) {
          plot <- plot + theme(axis.line.x = element_blank(),
                               axis.ticks.x = element_blank(),
                               axis.text.x = element_blank(),
                               axis.title.x = element_blank())
        }
      } else {
        if (myTheme == TRUE) {
          plot <- plot + theme(panel.border = element_rect(fill = NA, colour = "black",
                                                           linewidth = 0.8),
                               axis.title.x = element_blank(),
                               axis.ticks.x = element_blank(),
                               axis.ticks.y = element_blank(),
                               aspect.ratio = 1,
                               axis.title.y = element_blank(),
                               axis.line = element_blank(),
                               axis.text.x = element_blank(),
                               axis.text.y = element_blank(),
                               text = element_text(size = textsize))
        }
      }
      if (!blend) {
        plot <- plot + guides(color = NULL)
        cols.grad <- cols
        if (length(x = cols) == 1) {
          plot <- plot + scale_color_brewer(palette = cols)
        } else if (length(x = cols) > 1) {
          unique.feature.exp <- unique(data.plot[, feature])
          if (length(unique.feature.exp) == 1) {
            warn(message = paste0("All cells have the same value (",
                                  unique.feature.exp, ") of ", dQuote(x = feature)))
            if (unique.feature.exp == 0) { cols.grad <- cols[1] } else { cols.grad <- cols }
          }
          plot <- suppressMessages(expr = plot +
                                     scale_color_gradientn(colors = cols.grad,
                                                           guide = "colorbar"))
        }
      }
      if (!(is.null(x = keep.scale)) && keep.scale == "feature" && !blend) {
        max.feature.value <- max(data[, feature])
        min.feature.value <- min(data[, feature])
        plot <- suppressMessages(plot &
                                   scale_color_gradientn(colors = cols,
                                                         limits = c(min.feature.value,
                                                                    max.feature.value)))
      }
      if (coord.fixed) {
        plot <- plot + coord_fixed()
      }
      plots[[(length(x = features) * (i - 1)) + j]] <- plot
    }
  }
  if (isTRUE(x = blend)) {
    blend.legend <- BlendMap(color.matrix = color.matrix)
    for (ii in 1:length(x = levels(x = data$split))) {
      suppressMessages(expr = plots <- append(
        x = plots,
        values = list(blend.legend +
                        scale_y_continuous(sec.axis = dup_axis(name = ifelse(test = length(x = levels(x = data$split)) > 1,
                                                                              yes = levels(x = data$split)[ii],
                                                                              no = "")),
                                           expand = c(0, 0)) +
                        labs(x = features[1], y = features[2],
                             title = if (ii == 1) {
                               paste("Color threshold:", blend.threshold)
                             } else { NULL }) +
                        no.right),
        after = 4 * ii - 1))
    }
  }
  plots <- Filter(f = Negate(f = is.null), x = plots)
  if (is.null(x = ncol)) {
    ncol <- 2
    if (length(x = features) == 1) ncol <- 1
    if (length(x = features) > 6) ncol <- 3
    if (length(x = features) > 9) ncol <- 4
  }
  ncol <- ifelse(test = is.null(x = split.by) || isTRUE(x = blend),
                 yes = ncol, no = length(x = features))
  legend <- if (isTRUE(x = blend)) { "none" } else { split.by %iff% "none" }
  if (isTRUE(x = combine)) {
    if (by.col && !is.null(x = split.by) && !blend) {
      plots <- lapply(X = plots, FUN = function(x) {
        return(suppressMessages(expr = x + theme_cowplot() + ggtitle("") +
                                  scale_y_continuous(sec.axis = dup_axis(name = ""),
                                                     limits = ylims) + no.right))
      })
      nsplits <- length(x = levels(x = data$split))
      idx <- 1
      for (i in (length(x = features) * (nsplits - 1) + 1):(length(x = features) * nsplits)) {
        plots[[i]] <- suppressMessages(expr = plots[[i]] +
                                         scale_y_continuous(sec.axis = dup_axis(name = features[[idx]]),
                                                            limits = ylims) + no.right)
        idx <- idx + 1
      }
      idx <- 1
      for (i in which(x = 1:length(x = plots) %% length(x = features) == 1)) {
        plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
          theme(plot.title = element_text(hjust = 0.5))
        idx <- idx + 1
      }
      idx <- 1
      if (length(x = features) == 1) {
        for (i in 1:length(x = plots)) {
          plots[[i]] <- plots[[i]] + ggtitle(levels(x = data$split)[[idx]]) +
            theme(plot.title = element_text(hjust = 0.5))
          idx <- idx + 1
        }
        ncol <- 1
        nrow <- nsplits
      } else {
        nrow <- split.by %iff% length(x = levels(x = data$split))
      }
      plots <- plots[c(do.call(what = rbind,
                               args = split(x = 1:length(x = plots),
                                            f = ceiling(x = seq_along(along.with = 1:length(x = plots)) / length(x = features)))))]
      plots <- wrap_plots(plots, ncol = nrow, nrow = ncol)
      if (!is.null(x = legend) && legend == "none") {
        plots <- plots & NoLegend()
      }
    } else {
      plots <- wrap_plots(plots, ncol = ncol,
                          nrow = split.by %iff% length(x = levels(x = data$split)))
    }
    if (!is.null(x = legend) && legend == "none") {
      plots <- plots & NoLegend()
    }
    if (!(is.null(x = keep.scale)) && keep.scale == "all" && !blend) {
      max.feature.value <- max(data[, features])
      min.feature.value <- min(data[, features])
      plots <- suppressMessages(plots &
                                  scale_color_gradientn(colors = cols,
                                                        limits = c(min.feature.value,
                                                                   max.feature.value)))
    }
  }
  return(plots)
}

# Environment rebound to Seurat's namespace at load time (see R/zzz.R) so Seurat
# internals (SingleDimPlot, SetQuantile, CenterTitle, LabelClusters, BlendMatrix,
# BlendMap, BlendExpression, RandomName, NoLegend, IFeaturePlot, brewer.pal.info,
# %||%, %iff%, ...) resolve without :::
