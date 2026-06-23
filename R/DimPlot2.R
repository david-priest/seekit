#' @title DimPlot2
#' @description Seurat::DimPlot clone with a fixed publication theme:
#'   1:1 aspect ratio, panel border, no axis ticks/labels/lines, and
#'   adjustable legend point size. Identical signature to DimPlot for
#'   the parameters that matter; adds `legend_title`, `legendpointsize`,
#'   and `textsize`.
#'
#'   Uses Seurat's internal `SingleDimPlot`, `LabelClusters`, `FacetTheme`
#'   and `CenterTitle`. The function's environment is rebound to
#'   `asNamespace("Seurat")` at the end of this file so those internals
#'   resolve without explicit `:::`.
#'
#' @param object A Seurat object.
#' @param dims Two-length integer of reduction dims to plot. Default: c(1, 2)
#' @param cells Cells to include. Default: NULL (all)
#' @param cols Colour palette for `group.by` levels. Default: NULL
#' @param pt.size Point size. Default: NULL
#' @param reduction Name of the reduction to use. Default: NULL (Seurat default)
#' @param group.by Metadata column to colour by. Default: NULL (ident)
#' @param split.by Metadata column to facet by. Default: NULL
#' @param legend_title Legend title. Default: group.by
#' @param legendpointsize Override aes size in the legend. Default: 6
#' @param shape.by Metadata column to map to point shape. Default: NULL
#' @param order Order in which to plot cells. Default: NULL
#' @param shuffle Randomise plot order. Default: FALSE
#' @param seed RNG seed for shuffle. Default: 1
#' @param textsize Base text size. Default: 12
#' @param label Add cluster labels. Default: FALSE
#' @param label.size Label text size. Default: 4
#' @param label.color Label text colour. Default: "black"
#' @param label.box Draw box behind labels. Default: FALSE
#' @param repel Use ggrepel for labels. Default: FALSE
#' @param alpha Point alpha. Default: 1
#' @param cells.highlight Cells to highlight. Default: NULL
#' @param cols.highlight Highlight colour. Default: "#DE2D26"
#' @param sizes.highlight Highlight point size. Default: 1
#' @param na.value Colour for NA. Default: "grey50"
#' @param ncol Number of columns when combining. Default: NULL
#' @param combine Combine plots with patchwork::wrap_plots. Default: TRUE
#' @param raster Rasterise points. Default: NULL
#' @param raster.dpi Raster dpi as c(width, height). Default: c(512, 512)
#'
#' @return A ggplot (or patchwork object when multiple groups/splits).
#' @export
DimPlot2 <- function(object, dims = c(1, 2), cells = NULL, cols = NULL,
                     pt.size = NULL, reduction = NULL, group.by = NULL,
                     split.by = NULL, legend_title = group.by,
                     legendpointsize = 6, shape.by = NULL, order = NULL,
                     shuffle = FALSE, seed = 1, textsize = 12,
                     label = FALSE, label.size = 4, label.color = "black",
                     label.box = FALSE, repel = FALSE, alpha = 1,
                     cells.highlight = NULL, cols.highlight = "#DE2D26",
                     sizes.highlight = 1, na.value = "grey50",
                     ncol = NULL, combine = TRUE, raster = NULL,
                     raster.dpi = c(512, 512)) {
  if (!is_integerish(x = dims, n = 2L, finite = TRUE) || !all(dims > 0L)) {
    abort(message = "'dims' must be a two-length integer vector")
  }
  reduction <- reduction %||% DefaultDimReduc(object = object)
  cells <- cells %||% Cells(x = object, assay = DefaultAssay(object = object[[reduction]]))
  dims <- paste0(Key(object = object[[reduction]]), dims)
  orig.groups <- group.by
  group.by <- group.by %||% "ident"
  data <- FetchData(object = object, vars = c(dims, group.by),
                    cells = cells, clean = "project")
  group.by <- colnames(x = data)[3:ncol(x = data)]
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  if (!is.null(x = shape.by)) {
    data[, shape.by] <- object[[shape.by, drop = TRUE]]
  }
  if (!is.null(x = split.by)) {
    split <- FetchData(object = object, vars = split.by,
                       clean = TRUE)[split.by]
    data <- data[rownames(split), ]
    data[, split.by] <- split
  }
  if (isTRUE(x = shuffle)) {
    set.seed(seed = seed)
    data <- data[sample(x = 1:nrow(x = data)), ]
  }
  plots <- lapply(X = group.by, FUN = function(x) {
    plot <- SingleDimPlot(data = data[, c(dims, x, split.by, shape.by)],
                          dims = dims, col.by = x, cols = cols,
                          pt.size = pt.size, shape.by = shape.by, order = order,
                          alpha = alpha, label = FALSE,
                          cells.highlight = cells.highlight,
                          cols.highlight = cols.highlight,
                          sizes.highlight = sizes.highlight,
                          na.value = na.value, raster = raster,
                          raster.dpi = raster.dpi)
    if (label) {
      plot <- LabelClusters(plot = plot, id = x, repel = repel,
                            size = label.size, split.by = split.by,
                            box = label.box, color = label.color)
    }
    if (!is.null(x = split.by)) {
      plot <- plot + FacetTheme() +
        facet_wrap(facets = vars(!!sym(x = split.by)),
                   ncol = if (length(x = group.by) > 1 || is.null(x = ncol)) {
                     length(x = unique(x = data[, split.by]))
                   } else {
                     ncol
                   })
    }
    plot <- if (is.null(x = orig.groups)) {
      plot + labs(title = NULL)
    } else {
      plot + CenterTitle()
    }

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
                         text = element_text(size = textsize)) +
      guides(col = guide_legend(title = legend_title,
                                override.aes = list(alpha = 1,
                                                    size = legendpointsize)))
  })
  if (!is.null(x = split.by)) {
    ncol <- 1
  }
  if (combine) {
    plots <- wrap_plots(plots, ncol = orig.groups %iff% ncol)
  }
  return(plots)
}

# Environment rebound to Seurat's namespace at load time (see R/zzz.R) so Seurat
# internals (SingleDimPlot, LabelClusters, FacetTheme, CenterTitle, %iff%, etc.)
# resolve without :::.
