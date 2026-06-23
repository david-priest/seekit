#' @title Cairo-free rasterised points (drop-in for ggrastr::geom_point_rast)
#' @description A package-level `geom_point_rast` that rasterises via
#'   **scattermore** and NEVER touches Cairo. ggrastr's `geom_point_rast`
#'   uses a Cairo backend that fails ("Failed to create Cairo backend") or
#'   segfaults on large (hundreds-of-thousands-of-cell) UMAPs / scatterplots.
#'   This version draws a native rasterGrob with no offscreen device, so the
#'   full plot renders on any device (base `pdf()`, `ragg`, screen) without
#'   Cairo. It shadows `ggrastr::geom_point_rast` for all seekit
#'   plotters, so they get Cairo-free rasterisation automatically — no
#'   per-notebook shim required.
#'
#' @details ggrastr's `size` (mm) is translated to an approximate scattermore
#'   pixel radius (`pointsize`); the ggrastr-only `shape` / `raster.dpi` /
#'   `dev` arguments are accepted for call-compatibility and ignored.
#'
#' @param ... Passed to [scattermore::geom_scattermore()] (e.g. `position`).
#' @param size ggrastr point size (mm-ish); mapped to `pointsize` if the latter
#'   is not given. Default 0.4.
#' @param alpha Point alpha. Default 1.
#' @param shape,raster.dpi,dev ggrastr-only; accepted and ignored.
#' @param pointsize scattermore pixel radius. Default: `max(1.5, size * 3.5)`.
#' @param pixels scattermore raster size. Default `c(1200, 1200)`.
#' @return A ggplot2 layer.
#' @seealso [scattermore::geom_scattermore()]
#' @export
geom_point_rast <- function(..., size = 0.4, alpha = 1, shape = NULL,
                            raster.dpi = NULL, dev = NULL,
                            pointsize = NULL, pixels = c(1200, 1200)) {
  if (!requireNamespace("scattermore", quietly = TRUE))
    stop("scattermore is required for Cairo-free rasterisation. ",
         "install.packages('scattermore').")
  if (is.null(pointsize)) pointsize <- max(1.5, size * 3.5)
  scattermore::geom_scattermore(..., pointsize = pointsize, alpha = alpha,
                                pixels = pixels)
}
