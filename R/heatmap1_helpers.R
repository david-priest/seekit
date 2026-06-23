# heatmap1_helpers.R — CATALYST-free ports of the CustomFunctionsSept25 helpers
# that plotExprHeatmap1 (and plotExprHeatmapCol's sub=TRUE path) depend on:
# subSCE (already CATALYST-free), merging_ids, anno_clusters1.
# Bodies are David's own; only CATALYST's cluster_codes() -> .wl_cluster_codes().
# -----------------------------------------------------------------------------

# Per-sample downsample of an SCE (no CATALYST).
subSCE <- function(sce, n_cells = 1000, my_seed = 1234) {
  set.seed(my_seed)
  idx <- split(seq(ncol(sce)), sce$sample_id)
  idx <- lapply(idx, function(.) sample(., min(n_cells, length(.))))
  sce[, unlist(idx)]
}

# Build the cluster (+ optional merging) id/colour data.frame.
merging_ids <- function(x, k, m, k_pal, m_pal) {
  kids <- levels(x$cluster_id); nk <- length(kids)
  if (nk > length(k_pal)) k_pal <- grDevices::colorRampPalette(k_pal)(nk)
  k_pal <- k_pal[seq_len(nk)]; names(k_pal) <- kids
  df <- data.frame(cluster_id = kids); col <- list(cluster_id = k_pal)
  if (!is.null(m)) {
    i    <- match(kids, .wl_cluster_codes(x)[, k])
    mids <- droplevels(.wl_cluster_codes(x)[, m][i])
    nm   <- nlevels(mids)
    if (nm > length(m_pal)) m_pal <- grDevices::colorRampPalette(m_pal)(nm)
    m_pal <- m_pal[seq_len(nm)]; names(m_pal) <- levels(mids)
    df$merging_id <- mids; col$merging_id <- m_pal
  }
  dplyr::mutate_all(df, function(u) factor(u, unique(u)))
}

# Left ComplexHeatmap row-annotation for clusters (+ optional merging).
anno_clusters1 <- function(x, k, m, k_pal, m_pal, named_k_pal = NULL) {
  kids <- levels(x$cluster_id); nk <- length(kids)
  if (!is.null(named_k_pal)) {
    k_pal <- named_k_pal[kids]
  } else {
    if (nk > length(k_pal)) k_pal <- grDevices::colorRampPalette(k_pal)(nk)
    k_pal <- k_pal[seq_len(nk)]; names(k_pal) <- kids
  }
  df <- data.frame(cluster_id = kids); col <- list(cluster_id = k_pal)
  colout <<- col   # side-effect retained from the original (downstream reads it)
  if (!is.null(m)) {
    i    <- match(kids, .wl_cluster_codes(x)[, k])
    mids <- droplevels(.wl_cluster_codes(x)[, m][i])
    nm   <- nlevels(mids)
    if (nm > length(m_pal)) m_pal <- grDevices::colorRampPalette(m_pal)(nm)
    m_pal <- m_pal[seq_len(nm)]; names(m_pal) <- levels(mids)
    df$merging_id <- mids; col$merging_id <- m_pal
  }
  df <- dplyr::mutate_all(df, function(u) factor(u, unique(u)))
  ComplexHeatmap::rowAnnotation(df = df, col = col, gp = grid::gpar(col = "white"))
}
