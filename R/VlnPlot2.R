#' @title VlnPlot2
#' @description Seurat::VlnPlot drop-in with two additions on top of Seurat
#'   5.x's stock signature:
#'   - `alpha`: explicit alpha on the violin AND the jittered points
#'   - `raster`: when TRUE, jittered points are drawn with
#'     `ggrastr::geom_jitter_rast(raster.dpi = 600)` instead of
#'     `geom_jitter`. This keeps the violins as crisp vector geometry while
#'     rasterising the (typically thousands of) overlaid points so the PDF
#'     stays light and Illustrator-friendly.
#'
#'   `ExIPlot2` and `SingleExIPlot2` are dispatched copies of Seurat's
#'   internals carrying the same two additions; they're included in this
#'   file because `VlnPlot2` only works with the matched pair. Each of the
#'   three functions has its environment rebound to Seurat's namespace at
#'   the bottom of the file so the host of helpers they call (Layers,
#'   DefaultAssay, Cells, AutoPointSize, FetchData, MultiExIPlot,
#'   Interleave, InvertHex, Col2Hex, geom_split_violin, NoLegend,
#'   PackageCheck, hue_pal, wrap_plots, ...) resolve without explicit `:::`.
#'
#' @param object Seurat object.
#' @param features Features to plot (gene, metadata column, etc.).
#' @param cols Fill palette. Default: NULL
#' @param pt.size Point size for jittered points. 0 hides points. Default: NULL
#' @param alpha Alpha applied to BOTH violin fill and jittered points. Default: 1
#' @param idents Identities to include. Default: NULL
#' @param sort Sort identities by mean expression. Default: FALSE
#' @param assay Assay to pull data from. Default: NULL
#' @param group.by Metadata column to group cells by. Default: NULL
#' @param split.by Metadata column to split violins by. Default: NULL
#' @param adjust Violin bandwidth multiplier. Default: 1
#' @param y.max Manual y-axis max. Default: NULL
#' @param same.y.lims Use the same y-axis range across plots. Default: FALSE
#' @param log Log-scale the y-axis. Default: FALSE
#' @param ncol Patchwork columns. Default: NULL
#' @param slot Deprecated -- use layer.
#' @param layer Assay layer to pull from. Default: NULL (data)
#' @param split.plot Use single split-violin style for 2 groups. Default: FALSE
#' @param stack Stack feature plots vertically. Default: FALSE
#' @param combine Combine via patchwork. Default: TRUE
#' @param fill.by Which variable maps to fill. Default: "feature"
#' @param flip Flip x/y. Default: FALSE
#' @param add.noise Add tiny noise to avoid all-identical pile-ups. Default: TRUE
#' @param raster Rasterise jittered points via ggrastr::geom_jitter_rast.
#'   Default: NULL (auto: TRUE when >100k cells if ggrastr is installed).
#'
#' @return A ggplot (or patchwork) object.
#' @export
VlnPlot2 <- function(object, features, cols = NULL, pt.size = NULL, alpha = 1,
                     idents = NULL, sort = FALSE, assay = NULL, group.by = NULL,
                     split.by = NULL, adjust = 1, y.max = NULL,
                     same.y.lims = FALSE, log = FALSE, ncol = NULL,
                     slot = deprecated(), layer = NULL, split.plot = FALSE,
                     stack = FALSE, combine = TRUE, fill.by = "feature",
                     flip = FALSE, add.noise = TRUE, raster = NULL) {
  if (is_present(arg = slot)) {
    deprecate_soft(when = "5.0.0", what = "VlnPlot(slot = )",
                   with = "VlnPlot(layer = )")
    layer <- slot %||% layer
  }
  layer.set <- suppressWarnings(Layers(object = object,
                                       search = layer %||% "data"))
  if (is.null(layer) && length(layer.set) == 1 && layer.set == "scale.data") {
    warning("Default search for \"data\" layer yielded no results; utilizing \"scale.data\" layer instead.")
  }
  assay.name <- DefaultAssay(object)
  if (is.null(layer.set) & is.null(layer)) {
    warning("Default search for \"data\" layer in \"", assay.name,
            "\" assay yielded no results; utilizing \"counts\" layer instead.",
            call. = FALSE, immediate. = TRUE)
    layer.set <- Layers(object = object, search = "counts")
  }
  if (is.null(layer.set)) {
    stop("layer \"", layer, "\" is not found in assay: \"", assay.name, "\"")
  } else {
    layer <- layer.set
  }
  if (!is.null(x = split.by) &
      getOption(x = "Seurat.warn.vlnplot.split", default = TRUE)) {
    message("The default behaviour of split.by has changed.\n",
            "Separate violin plots are now plotted side-by-side.\n",
            "To restore the old behaviour of a single split violin,\n",
            "set split.plot = TRUE.\n      \nThis message will be shown once per session.")
    options(Seurat.warn.vlnplot.split = FALSE)
  }
  return(ExIPlot2(object = object,
                  type = ifelse(test = split.plot, yes = "splitViolin", no = "violin"),
                  features = features, idents = idents, ncol = ncol, sort = sort,
                  assay = assay, y.max = y.max, same.y.lims = same.y.lims,
                  adjust = adjust, pt.size = pt.size, alpha = alpha, cols = cols,
                  group.by = group.by, split.by = split.by, log = log,
                  layer = layer, stack = stack, combine = combine,
                  fill.by = fill.by, flip = flip, add.noise = add.noise,
                  raster = raster))
}

# ── ExIPlot2 ─────────────────────────────────────────────────────────────────
ExIPlot2 <- function(object, features, type = "violin", idents = NULL,
                     ncol = NULL, sort = FALSE, assay = NULL, y.max = NULL,
                     same.y.lims = FALSE, adjust = 1, cols = NULL, pt.size = 0,
                     alpha = 1, group.by = NULL, split.by = NULL, log = FALSE,
                     slot = deprecated(), layer = "data", stack = FALSE,
                     combine = TRUE, fill.by = NULL, flip = FALSE,
                     add.noise = TRUE, raster = NULL) {
  if (is_present(arg = slot)) {
    layer <- layer %||% slot
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  cells <- Cells(x = object, assay = NULL)
  if (isTRUE(x = stack)) {
    if (!is.null(x = ncol)) {
      warning("'ncol' is ignored with 'stack' is TRUE", call. = FALSE, immediate. = TRUE)
    }
    if (!is.null(x = y.max)) {
      warning("'y.max' is ignored when 'stack' is TRUE", call. = FALSE, immediate. = TRUE)
    }
  } else {
    ncol <- ncol %||% ifelse(test = length(x = features) > 9, yes = 4,
                             no = min(length(x = features), 3))
  }
  if (!is.null(x = idents)) {
    cells <- intersect(x = names(x = Idents(object = object)[Idents(object = object) %in% idents]),
                       y = cells)
  }
  data <- FetchData(object = object, vars = features, layer = layer, cells = cells)
  pt.size <- pt.size %||% AutoPointSize(data = object)
  features <- colnames(x = data)
  data <- data[cells, , drop = FALSE]
  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  } else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- factor(x = idents)
  }
  if (is.null(x = split.by)) {
    split <- NULL
  } else {
    split <- FetchData(object, split.by)[cells, split.by]
    if (!is.factor(x = split)) {
      split <- factor(x = split)
    }
    if (is.null(x = cols)) {
      cols <- hue_pal()(length(x = levels(x = idents)))
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    } else if (length(x = cols) == 1 && cols == "interaction") {
      split <- interaction(idents, split)
      cols <- hue_pal()(length(x = levels(x = idents)))
    } else {
      cols <- Col2Hex(cols)
    }
    if (length(x = cols) < length(x = levels(x = split))) {
      cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    }
    cols <- rep_len(x = cols, length.out = length(x = levels(x = split)))
    names(x = cols) <- levels(x = split)
    if ((length(x = cols) > 2) & (type == "splitViolin")) {
      warning("Split violin is only supported for <3 groups, using multi-violin.")
      type <- "violin"
    }
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data)
  }
  if (isTRUE(x = stack)) {
    return(MultiExIPlot(type = type, data = data, idents = idents, split = split,
                        sort = sort, same.y.lims = same.y.lims, adjust = adjust,
                        cols = cols, pt.size = pt.size, log = log, fill.by = fill.by,
                        add.noise = add.noise, flip = flip))
  }
  plots <- lapply(X = features, FUN = function(x) {
    return(SingleExIPlot2(type = type, data = data[, x, drop = FALSE],
                          idents = idents, split = split, sort = sort,
                          y.max = y.max, adjust = adjust, cols = cols,
                          pt.size = pt.size, alpha = alpha, log = log,
                          add.noise = add.noise, raster = raster))
  })
  label.fxn <- switch(EXPR = type,
                      violin = if (stack) { xlab } else { ylab },
                      splitViolin = if (stack) { xlab } else { ylab },
                      ridge = xlab,
                      stop("Unknown ExIPlot type ", type, call. = FALSE))
  for (i in 1:length(x = plots)) {
    key <- paste0(unlist(x = strsplit(x = features[i], split = "_"))[1], "_")
    obj <- names(x = which(x = Key(object = object) == key))
    if (length(x = obj) == 1) {
      if (inherits(x = object[[obj]], what = "DimReduc")) {
        plots[[i]] <- plots[[i]] + label.fxn(label = "Embeddings Value")
      } else if (inherits(x = object[[obj]], what = "Assay") ||
                 inherits(x = object[[obj]], what = "Assay5")) {
        next
      } else {
        warning("Unknown object type ", class(x = object), immediate. = TRUE, call. = FALSE)
        plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
      }
    } else if (!features[i] %in% rownames(x = object)) {
      plots[[i]] <- plots[[i]] + label.fxn(label = NULL)
    }
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = ncol)
    if (length(x = features) > 1) {
      plots <- plots & NoLegend()
    }
  }
  return(plots)
}

# ── SingleExIPlot2 ───────────────────────────────────────────────────────────
SingleExIPlot2 <- function(data, idents, split = NULL, type = "violin",
                           sort = FALSE, y.max = NULL, adjust = 1, pt.size = 0,
                           alpha = 1, cols = NULL, seed.use = 42, log = FALSE,
                           add.noise = TRUE, raster = NULL) {
  if (!is.null(x = raster) && isTRUE(x = raster)) {
    if (!requireNamespace("scattermore", quietly = TRUE)) {
      stop("Please install scattermore from CRAN to enable (Cairo-free) rasterization.")
    }
  }
  if (requireNamespace("scattermore", quietly = TRUE)) {
    if ((nrow(x = data) > 1e+05) & is.null(x = raster)) {
      message("Rasterizing points since number of points exceeds 100,000.",
              "\nTo disable this behavior set `raster=FALSE`")
      raster <- TRUE
    }
  }
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  if (!is.data.frame(x = data) || ncol(x = data) != 1) {
    stop("'SingleExIPlot requires a data frame with 1 column")
  }
  feature <- colnames(x = data)
  data$ident <- idents
  if ((is.character(x = sort) && nchar(x = sort) > 0) || sort) {
    data$ident <- factor(x = data$ident,
                         levels = names(x = rev(x = sort(x = tapply(X = data[, feature],
                                                                    INDEX = data$ident,
                                                                    FUN = mean),
                                                         decreasing = grepl(pattern = paste0("^", tolower(x = sort)),
                                                                            x = "decreasing")))))
  }
  if (log) {
    noise <- rnorm(n = length(x = data[, feature])) / 200
    data[, feature] <- data[, feature] + 1
  } else {
    noise <- rnorm(n = length(x = data[, feature])) / 1e+05
  }
  if (!add.noise) {
    noise <- noise * 0
  }
  if (all(data[, feature] == data[, feature][1])) {
    warning(paste0("All cells have the same value of ", feature, "."))
  } else {
    data[, feature] <- data[, feature] + noise
  }
  axis.label <- "Expression Level"
  y.max <- y.max %||% max(data[, feature][is.finite(x = data[, feature])])
  if (type == "violin" && !is.null(x = split)) {
    data$split <- split
    vln.geom <- geom_violin
    fill <- "split"
  } else if (type == "splitViolin" && !is.null(x = split)) {
    data$split <- split
    vln.geom <- geom_split_violin
    fill <- "split"
    type <- "violin"
  } else {
    vln.geom <- geom_violin
    fill <- "ident"
  }
  switch(EXPR = type,
         violin = {
           x <- "ident"
           y <- paste0("`", feature, "`")
           xlab <- "Identity"
           ylab <- axis.label
           geom <- list(vln.geom(scale = "width", adjust = adjust, trim = TRUE),
                        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
           if (is.null(x = split)) {
             jitter <- if (isTRUE(x = raster)) {
               # Cairo-free rasterised jitter via scattermore (NOT ggrastr/Cairo)
               scattermore::geom_scattermore(position = position_jitter(height = 0),
                                             pointsize = max(1.5, pt.size * 3.5),
                                             alpha = alpha, show.legend = FALSE)
             } else {
               geom_jitter(height = 0, size = pt.size, alpha = alpha,
                           show.legend = FALSE)
             }
           } else {
             jitter <- if (isTRUE(x = raster)) {
               # Cairo-free rasterised jitter via scattermore (NOT ggrastr/Cairo)
               scattermore::geom_scattermore(position = position_jitterdodge(jitter.width = 0.4,
                                                                             dodge.width = 0.9),
                                             pointsize = max(1.5, pt.size * 3.5),
                                             alpha = alpha, show.legend = FALSE)
             } else {
               geom_jitter(position = position_jitterdodge(jitter.width = 0.4,
                                                           dodge.width = 0.9),
                           size = pt.size, alpha = alpha, show.legend = FALSE)
             }
           }
           log.scale <- scale_y_log10()
           axis.scale <- ylim
         },
         ridge = {
           x <- paste0("`", feature, "`")
           y <- "ident"
           xlab <- axis.label
           ylab <- "Identity"
           geom <- list(geom_density_ridges(scale = 4), theme_ridges(),
                        scale_y_discrete(expand = c(0.01, 0)),
                        scale_x_continuous(expand = c(0, 0)))
           jitter <- geom_jitter(width = 0, size = pt.size, alpha = alpha,
                                 show.legend = FALSE)
           log.scale <- scale_x_log10()
           axis.scale <- function(...) { invisible(x = NULL) }
         },
         stop("Unknown plot type: ", type))
  plot <- ggplot(data = data,
                 mapping = aes_string(x = x, y = y, fill = fill)[c(2, 3, 1)]) +
    labs(x = xlab, y = ylab, title = feature, fill = NULL) +
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5))
  plot <- do.call(what = "+", args = list(plot, geom))
  plot <- plot + if (log) {
    log.scale
  } else {
    axis.scale(min(data[, feature]), y.max)
  }
  if (pt.size > 0) {
    plot <- plot + jitter
  }
  if (!is.null(x = cols)) {
    if (!is.null(x = split)) {
      idents <- unique(x = as.vector(x = data$ident))
      splits <- unique(x = as.vector(x = data$split))
      labels <- if (length(x = splits) == 2) {
        splits
      } else {
        unlist(x = lapply(X = idents, FUN = function(pattern, x) {
          x.mod <- gsub(pattern = paste0(pattern, "."),
                        replacement = paste0(pattern, ": "), x = x, fixed = TRUE)
          x.keep <- grep(pattern = ": ", x = x.mod, fixed = TRUE)
          x.return <- x.mod[x.keep]
          names(x = x.return) <- x[x.keep]
          return(x.return)
        }, x = unique(x = as.vector(x = data$split))))
      }
      if (is.null(x = names(x = labels))) {
        names(x = labels) <- labels
      }
    } else {
      labels <- levels(x = droplevels(data$ident))
    }
    plot <- plot + scale_fill_manual(values = cols, labels = labels)
  }
  return(plot)
}

# Environments rebound to Seurat's namespace at load time (see R/zzz.R) so each
# function can reach Seurat internals (Layers, AutoPointSize, MultiExIPlot,
# Interleave, InvertHex, Col2Hex, geom_split_violin, NoLegend, PackageCheck,
# hue_pal, %||%, ...) without :::
