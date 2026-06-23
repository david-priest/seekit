# internals.R
# -----------------------------------------------------------------------------
# Clean-room MIT reimplementations of the CATALYST internals/accessors that the
# winglab CyTOF plotters depend on. Written from the documented behaviour and
# the golden-reference I/O (see capture_golden_reference.R) — NOT from CATALYST
# source. No CATALYST dependency. Independent expression.
#
# Functional contracts only (these are the "method of operation" — aggregate by
# group, quantile-rescale, resolve a marker class — which are not protectable);
# the implementations below are original.
# -----------------------------------------------------------------------------
# Deps (already in seekit): SingleCellExperiment, SummarizedExperiment,
# S4Vectors, matrixStats.

# ---- accessors --------------------------------------------------------------

# Per-CATALYST-convention an SCE carries: metadata(x)$cluster_codes (a data.frame,
# one row per base SOM cluster, one column per clustering resolution), a per-cell
# base clustering in colData(x)$cluster_id, sample labels in colData(x)$sample_id,
# and the experiment-info table in metadata(x)$experiment_info.

.wl_cluster_codes <- function(x) S4Vectors::metadata(x)$cluster_codes

.wl_sample_ids <- function(x) x$sample_id

.wl_ei <- function(x) S4Vectors::metadata(x)$experiment_info

# cluster_ids(x, k): map each cell's base cluster to its label at resolution k.
.wl_cluster_ids <- function(x, k = NULL) {
  if (is.null(k)) return(x$cluster_id)
  k  <- .wl_check_k(x, k)
  cc <- .wl_cluster_codes(x)
  factor(cc[[k]][as.numeric(x$cluster_id)], levels = levels(cc[[k]]))
}

# [#7] Resolve per-cell cluster labels for `k` the SAME way every plotter does:
# prefer a colData column named `k` (what plotExprHeatmap*/plotDR2 use when k is
# in colData), else fall back to the cluster_codes mapping. Keeping palette and
# plotters on one resolver prevents the colData-vs-cluster_codes drift that
# silently dropped annotation colours (the missing-circles bug).
.wl_resolve_cluster_ids <- function(x, k) {
  cd <- SummarizedExperiment::colData(x)
  if (!is.null(k) && length(k) == 1 && k %in% colnames(cd)) factor(cd[[k]])
  else .wl_cluster_ids(x, k)
}

# ---- validators -------------------------------------------------------------

# These mirror CATALYST's .check_* contracts: NULL-tolerant and returning TRUE so
# they compose inside stopifnot(); .wl_check_k returns the resolved clustering name.

.wl_check_sce <- function(x, needs_clustering = FALSE) {
  stopifnot(methods::is(x, "SingleCellExperiment"), !is.null(x$sample_id))
  if (needs_clustering)
    stopifnot(!is.null(x$cluster_id), !is.null(.wl_cluster_codes(x)))
  TRUE
}

# resolution name ("meta20"/"merging1") or numeric k -> "meta{k}"; NULL -> first.
.wl_check_k <- function(x, k = NULL) {
  cc <- .wl_cluster_codes(x)
  if (is.null(k)) return(colnames(cc)[1])
  if (is.numeric(k)) k <- paste0("meta", k)
  k <- as.character(k)
  if (length(k) != 1 || !k %in% colnames(cc))
    stop("'k' = '", paste(k, collapse = ","), "' is not a clustering. Available: ",
         paste(colnames(cc), collapse = ", "))
  k
}

.wl_check_assay <- function(x, y) {
  stopifnot(length(y) == 1, is.character(y),
            sum(y == SummarizedExperiment::assayNames(x)) == 1)
  TRUE
}

.wl_check_pal <- function(x, n = 2) {
  if (is.null(x)) return(TRUE)
  stopifnot(is.character(x), length(x) >= n)
  if (is.null(tryCatch(grDevices::col2rgb(x), error = function(e) NULL)))
    stop("Invalid colour palette.")
  TRUE
}

.wl_check_cd_factor <- function(x, y, n = 1) {
  if (is.null(y)) return(TRUE)
  if (!is.null(n)) stopifnot(length(y) == n)
  cd <- SummarizedExperiment::colData(x)
  stopifnot(is.character(y), all(y %in% names(cd)),
            !vapply(cd[y], is.numeric, logical(1)))
  TRUE
}

# ---- feature resolution -----------------------------------------------------

# fs: "type"/"state"/"none" -> by marker_class; otherwise an explicit gene vector.
.wl_get_features <- function(x, fs = NULL) {
  if (is.null(fs)) fs <- "type"
  if (length(fs) == 1 && fs %in% c("type", "state", "none")) {
    if (fs == "none") return(rownames(x))
    mc <- SummarizedExperiment::rowData(x)$marker_class
    return(rownames(x)[!is.na(mc) & mc == fs])
  }
  if (!all(fs %in% rownames(x)))
    stop("Unknown features: ", paste(setdiff(fs, rownames(x)), collapse = ", "))
  fs
}

# ---- aggregation ------------------------------------------------------------

# Aggregate an assay (markers x cells) over one or two grouping columns.
#   by length 1 -> matrix (markers x groups)
#   by length 2 -> named list, one per level of by[1], each a markers x
#                  levels(by[2]) matrix (matches CATALYST's 2D aggregation shape)
.wl_agg <- function(x, by = c("cluster_id", "sample_id"),
                    fun = c("median", "mean", "sum"), assay = "exprs") {
  fun <- match.arg(fun)
  by  <- match.arg(by, several.ok = TRUE)
  y   <- as.matrix(SummarizedExperiment::assay(x, assay))   # markers x cells
  rowFun <- switch(fun,
    median = function(m) matrixStats::rowMedians(m, useNames = FALSE),
    mean   = function(m) base::rowMeans(m),
    sum    = function(m) base::rowSums(m))

  if (length(by) == 1) {
    grp <- factor(SummarizedExperiment::colData(x)[[by]])
    idx <- split(seq_len(ncol(y)), grp)
    res <- vapply(idx, function(j) rowFun(y[, j, drop = FALSE]), numeric(nrow(y)))
    rownames(res) <- rownames(y)
    return(res)
  }

  g1   <- factor(SummarizedExperiment::colData(x)[[by[1]]])
  g2   <- factor(SummarizedExperiment::colData(x)[[by[2]]])
  lev2 <- levels(g2)
  cells1 <- split(seq_len(ncol(y)), g1)              # cells per by[1] level
  lapply(cells1, function(j) {
    sub  <- y[, j, drop = FALSE]
    idx2 <- split(seq_along(j), factor(g2[j], levels = lev2))   # keep all by[2] levels
    m <- vapply(idx2, function(jj) {
      if (!length(jj)) rep(NA_real_, nrow(y)) else rowFun(sub[, jj, drop = FALSE])
    }, numeric(nrow(y)))
    rownames(m) <- rownames(y)
    m
  })
}

# ---- expression scaling -----------------------------------------------------

# Quantile-clip to [q, 1-q] along `margin` (1 = rows, 2 = cols) then rescale to
# [0, 1]. Degenerate ranges (hi == lo) map to 0.
.wl_scale_exprs <- function(x, margin = 1, q = 0.01) {
  x <- as.matrix(x)
  qfun <- function(v) stats::quantile(v, probs = c(q, 1 - q), na.rm = TRUE)
  if (margin == 1) {
    rng <- apply(x, 1, qfun); lo <- rng[1, ]; hi <- rng[2, ]
    sc  <- (x - lo) / (hi - lo)            # lo/hi recycle down each column (per-row)
  } else {
    rng <- apply(x, 2, qfun); lo <- rng[1, ]; hi <- rng[2, ]
    sc  <- sweep(sweep(x, 2, lo, "-"), 2, hi - lo, "/")
  }
  sc[sc < 0] <- 0; sc[sc > 1] <- 1
  sc[!is.finite(sc)] <- 0
  sc
}

# ---- heatmap factor annotations ---------------------------------------------

# Build a ComplexHeatmap annotation from the colData factors that are constant
# within each sample (i.e. sample-level metadata such as condition/patient).
# `which` may be a character vector restricting which factor(s) to show.
# Returns a HeatmapAnnotation, or NULL if no eligible factor.
.wl_anno_factors <- function(x, ids, which, type = c("row", "column")) {
  type <- match.arg(type)
  cd <- as.data.frame(SummarizedExperiment::colData(x), check.names = FALSE)
  cd <- cd[, !vapply(cd, is.numeric, logical(1)), drop = FALSE]
  cd[] <- lapply(cd, function(v) droplevels(factor(v)))
  # columns with exactly one level within every sample -> sample-level metadata
  one_per_sample <- vapply(cd, function(v)
    all(tapply(v, cd$sample_id, function(s) nlevels(droplevels(factor(s)))) == 1),
    logical(1))
  keep <- setdiff(names(cd)[one_per_sample], c("sample_id", "cluster_id"))
  if (is.character(which)) keep <- intersect(keep, which)
  if (length(keep) == 0) return(NULL)

  m  <- match(ids, cd$sample_id)
  df <- cd[m, keep, drop = FALSE]
  lvls  <- lapply(df, levels)
  nlvls <- vapply(lvls, length, numeric(1))
  pal <- RColorBrewer::brewer.pal(8, "Set3")[-2]
  if (any(nlvls > length(pal))) pal <- grDevices::colorRampPalette(pal)(max(nlvls))
  cols <- stats::setNames(lapply(names(df), function(i) {
    u <- pal[seq_len(nlvls[i])]; names(u) <- lvls[[i]]; u
  }), names(df))
  ComplexHeatmap::HeatmapAnnotation(which = type, df = df, col = cols,
                                    gp = grid::gpar(col = "white"))
}

# ---- point shapes -----------------------------------------------------------

# Map the levels of a `shape_by` factor to point shapes (preferred solid/line
# shapes first, then fill out from the remaining pch values). NULL if no
# shape_by, or if >17 levels (too many to distinguish).
.wl_get_shapes <- function(x, shape_by) {
  if (is.null(shape_by)) return(NULL)
  n    <- nlevels(x[[shape_by]])
  pref <- c(16, 17, 15, 3, 7, 8)
  if (n > 18) {
    message("At most 17 shapes are supported but ", n,
            " are required; setting 'shape_by' to NULL.")
    return(NULL)
  }
  if (n > length(pref)) pref <- c(pref, setdiff(c(0:15, 18), pref))
  pref[seq_len(n)]
}

# ---- Cairo-free rasterised points -------------------------------------------

# Drop-in replacement for ggrastr::geom_point_rast that NEVER touches Cairo.
# ggrastr's Cairo backend fails ("Failed to create Cairo backend") / segfaults on
# large (~370k-cell) UMAPs; scattermore rasterises natively with no Cairo and no
# offscreen device. Translates ggrastr's `size` to a scattermore pixel radius;
# the ggrastr-only `shape`/`raster.dpi`/`dev` args are accepted and ignored.
.wl_geom_point_rast <- function(..., size = 0.4, alpha = 1, shape = NULL,
                                raster.dpi = NULL, dev = NULL,
                                pointsize = NULL, pixels = c(1200, 1200)) {
  if (is.null(pointsize)) pointsize <- max(1.5, size * 3.5)
  scattermore::geom_scattermore(..., pointsize = pointsize, alpha = alpha, pixels = pixels)
}

# ---- cluster palette --------------------------------------------------------

# Conventional 30-colour qualitative palette (ColorBrewer Paired/Set1-derived
# hues + extensions) used for cluster colouring. A list of colours, freely
# replaceable — change here to restyle every plot.
.wl_cluster_cols <- c(
  "#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", "#B17BA6",
  "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", "#33A02C", "#B2DF8A",
  "#55A1B1", "#8DD3C7", "#A6761D", "#E6AB02", "#7570B3", "#BEAED4",
  "#666666", "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3",
  "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
