#' @title Compensate a SingleCellExperiment (CATALYST-free), with a progress bar
#' @description Clean-room, CATALYST-free re-implementation of `CATALYST::compCytof()`
#'   for mass-cytometry compensation. Applies a spillover matrix `sm` to an assay
#'   (default `counts`), writes the compensated values back (overwriting the assay, or
#'   into `compcounts`), and optionally adds the arcsinh `exprs`. The per-cell NNLS solve
#'   — the slow part on big (e.g. un-debarcoded, concatenated) files — shows a
#'   Seurat-style progress bar via `pbapply`.
#'
#'   Differences vs CATALYST: the spillover matrix is aligned to the data channels by
#'   **name only** (no isotope-impurity inference). This is correct when `sm` already
#'   covers the channels you compensate — e.g. a `computeSpillmat` / `update_spill_matrix2`
#'   matrix over the panel + barcode channels. Channels present in the data but not in
#'   `sm` are left uncompensated (identity diagonal); channels in `sm` but not the data
#'   are dropped. Because `sm` is expanded to the full data-channel set (identity for
#'   unknowns), there is no dimension mismatch even when the SCE still carries non-mass
#'   channels (Time, Gaussian) — they pass through untouched.
#'
#'   `cofactor` is treated as a single scalar (the common case); per-channel cofactors
#'   are not reproduced.
#' @param x A `SingleCellExperiment` (rows = channels; `rowData(x)$channel_name` set, as
#'   produced by `prepData2` / `CATALYST::prepData`).
#' @param sm Spillover matrix (rows = emitting, cols = receiving channels), named by
#'   channel on both dims. If `NULL`, taken from `metadata(x)$spillover_matrix`.
#' @param method `"nnls"` (per-cell non-negative least squares; default) or `"flow"`
#'   (vectorised linear unmixing via `flowCore::compensate`).
#' @param assay Assay to compensate, Default: 'counts'.
#' @param overwrite Overwrite `assay` (and `exprs`) in place; else write to
#'   `compcounts` / `compexprs`, Default: TRUE.
#' @param transform Add an arcsinh assay computed from the compensated counts,
#'   Default: TRUE.
#' @param cofactor arcsinh cofactor; if `NULL`, taken from `int_metadata(x)$cofactor`,
#'   Default: NULL.
#' @param verbose Show a progress bar on the NNLS solve, Default: TRUE.
#' @param cl Passed to `pbapply::pbapply` for parallelism (integer cores or a cluster);
#'   `NULL` = serial. Only used when `verbose = TRUE`, Default: NULL.
#' @return The `SingleCellExperiment` with compensated assay(s); cell order preserved.
#' @export
#' @importFrom SummarizedExperiment assay assay<- assayNames rowData
#' @importFrom SingleCellExperiment int_metadata int_metadata<-
#' @importFrom S4Vectors metadata
#' @importFrom flowCore flowFrame compensate exprs
#' @importFrom nnls nnls
#' @importFrom pbapply pbapply
#' @importFrom methods is
compCytof2 <- function(x, sm = NULL, method = c("nnls", "flow"), assay = "counts",
                       overwrite = TRUE, transform = TRUE, cofactor = NULL,
                       verbose = TRUE, cl = NULL) {
  method <- match.arg(method)
  stopifnot(methods::is(x, "SingleCellExperiment"),
            assay %in% SummarizedExperiment::assayNames(x),
            is.logical(overwrite), length(overwrite) == 1,
            is.logical(transform),  length(transform)  == 1)

  if (is.null(sm)) sm <- S4Vectors::metadata(x)$spillover_matrix
  if (is.null(sm)) stop("compCytof2: no `sm`, and none in metadata(x)$spillover_matrix.", call. = FALSE)
  if (!is.matrix(sm)) sm <- as.matrix(sm)
  if (is.null(cofactor)) cofactor <- SingleCellExperiment::int_metadata(x)$cofactor
  if (transform && (is.null(cofactor) || !is.numeric(cofactor)))
    stop("compCytof2: `transform = TRUE` needs a numeric `cofactor` (arg or int_metadata).", call. = FALSE)

  # work in channel space (sm is channel-named); restore marker rownames at the end
  chs0 <- rownames(x)
  chs  <- as.character(SummarizedExperiment::rowData(x)$channel_name)
  if (any(is.na(chs))) chs[is.na(chs)] <- chs0[is.na(chs)]
  rownames(x) <- chs

  sm <- .wl_adapt_sm(sm, chs)                       # expand/contract sm to the data channels

  y <- switch(method,
    flow = {
      a  <- as.matrix(SummarizedExperiment::assay(x, assay))
      ff <- flowCore::flowFrame(t(a))
      t(flowCore::exprs(flowCore::compensate(ff, sm)))
    },
    nnls = {
      smt <- t(sm)
      a   <- as.matrix(SummarizedExperiment::assay(x, assay))
      fn  <- function(u) nnls::nnls(smt, u)$x
      out <- if (verbose) pbapply::pbapply(a, 2, fn, cl = cl) else apply(a, 2, fn)
      rownames(out) <- rownames(a)
      out
    })

  aout_c <- if (overwrite) assay else "compcounts"
  SummarizedExperiment::assay(x, aout_c, withDimnames = FALSE) <- y
  if (transform) {
    aout_e <- if (overwrite) "exprs" else "compexprs"
    SummarizedExperiment::assay(x, aout_e, withDimnames = FALSE) <-
      asinh(SummarizedExperiment::assay(x, aout_c) / cofactor[1])
    SingleCellExperiment::int_metadata(x)$cofactor <- cofactor
  }
  rownames(x) <- chs0
  x
}

# Name-based spillover alignment to a target channel set: CATALYST's adaptSpillmat()
# minus the isotope-impurity inference. Identity diagonal for channels with no estimated
# spillover; values copied where channels overlap `sm`.
.wl_adapt_sm <- function(sm, chs) {
  out <- diag(length(chs)); dimnames(out) <- list(chs, chs)
  r  <- intersect(rownames(sm), chs)
  cc <- intersect(colnames(sm), chs)
  if (length(r) && length(cc)) out[r, cc] <- sm[r, cc]
  diag(out) <- 1
  out
}
