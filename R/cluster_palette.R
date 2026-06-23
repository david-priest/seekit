# cluster_palette.R — CATALYST-free rewrite of cluster_palette.
# Body is David's own; only CATALYST::cluster_ids -> .wl_cluster_ids and
# CATALYST:::.cluster_cols -> .wl_cluster_cols.
# -----------------------------------------------------------------------------
cluster_palette <- function(sce, k = "merging1", palette_fn = NULL) {
  # [#7] resolve clusters the same way the plotters do (colData-preferred) so the
  # palette always covers every cluster the heatmap/UMAP will draw.
  cluster_levels <- levels(.wl_resolve_cluster_ids(sce, k))
  n <- length(cluster_levels)
  if (n == 0) stop("cluster_palette: no levels found for k = '", k, "'")

  if (is.null(palette_fn)) {
    pal <- .wl_cluster_cols
    pal <- if (n > length(pal)) rep_len(pal, n) else pal[seq_len(n)]
  } else {
    pal <- palette_fn(n)
    if (length(pal) != n)
      stop("cluster_palette: palette_fn returned ", length(pal),
           " colours but ", n, " were needed.")
  }
  stats::setNames(pal, cluster_levels)
}
