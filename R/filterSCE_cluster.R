# filterSCE_cluster
#
# CATALYST-free filter of a SingleCellExperiment by cluster membership at a given
# clustering resolution `k`. This is the MIT replacement for the CATALYST-derived
# pattern `filterSCEB(x, cluster_id %in% ..., k = metalevel)`: it resolves the
# per-cell cluster label via the package internals (.wl_resolve_cluster_ids), so
# it works whether `k` is a colData column or a cluster_codes resolution such as
# "meta20" / "merging2". Sibling of filterSCE_simple2 (which filters by a plain
# colData column). 2026-06 seekit migration of the CMV CyTOF pipeline.
#
# @param sce      A SingleCellExperiment.
# @param clusters Cluster label(s) to keep (or drop, if exclude = TRUE).
# @param k        Clustering resolution: a colData column name, a cluster_codes
#                 column ("meta20", "merging1", ...), or numeric (-> "meta{k}").
# @param exclude  If TRUE, keep cells NOT in `clusters`. Default FALSE.

filterSCE_cluster <- function(sce, clusters, k, exclude = FALSE) {
  stopifnot(methods::is(sce, "SingleCellExperiment"))
  ids       <- as.character(.wl_resolve_cluster_ids(sce, k))
  condition <- ids %in% as.character(clusters)
  if (exclude) condition <- !condition
  condition[is.na(condition)] <- FALSE

  sce_filtered <- sce[, condition]

  # Drop unused levels in all factor columns in the colData (matches filterSCE_simple2).
  is_factor <- sapply(colData(sce_filtered), is.factor)
  colData(sce_filtered)[is_factor] <- lapply(colData(sce_filtered)[is_factor], droplevels)

  return(sce_filtered)
}
