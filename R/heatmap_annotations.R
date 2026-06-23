# heatmap_annotations
#
# Internal ComplexHeatmap row-annotation helpers used by plotExprHeatmap1() and
# plotExprHeatmapCol(): coloured cluster dots on the left (anno_clusters2 ->
# circle_anno) and the per-cluster count / percentage annotation on the right
# (anno_counts1). Extracted from "CustomFunctions not Annotated New 250131.Rmd"
# — the same source the heatmap functions were migrated from — which had left
# these helpers behind, so both plotters errored with
# `could not find function "anno_clusters2"`.
#
# These rely on ComplexHeatmap / grid being attached (the analysis notebooks
# library(ComplexHeatmap)); the calling plotters make the same assumption.

# Left annotation: one filled colour dot per heatmap row, coloured by cluster.
anno_clusters2 <- function(cluster_order, k_pal) {
  ComplexHeatmap::rowAnnotation(
    annotation_function = function(index) {
      circle_anno(index, k_pal, cluster_order)
    },
    show_annotation_name = FALSE,
    gp = gpar(col = "white")
  )
}

# Draw the coloured dots for anno_clusters2 (called per heatmap render).
circle_anno <- function(index, k_pal, cluster_order) {
  n <- length(index)
  # Map each row to its cluster's colour via the row order.
  colors <- k_pal[cluster_order[index]]
  grid.points(x = rep(0.5, n),
              y = unit((n:1 - 0.5) / n, "npc"),
              pch = 16,
              size = unit(8, "mm"),
              gp = gpar(col = colors))
}

# Right annotation: per-cluster cell-count barplot, optionally labelled with the
# percentage of total cells in each cluster.
anno_counts1 <- function(x, perc) {
  ns <- table(x)
  fq <- round(ns / sum(ns) * 100, 2)
  if (perc) {
    txt <- paste(fq, "%")
    foo <- row_anno_text(txt, just = "center", gp = gpar(fontsize = 12),
                         location = unit(0.5, "npc"))
  } else {
    foo <- NULL
  }
  ComplexHeatmap::rowAnnotation(
    n_cells = row_anno_barplot(x = as.matrix(ns), width = unit(2, "cm"),
                               gp = gpar(fill = "grey", col = "white"),
                               border = FALSE, axis = TRUE, bar_width = 0.8),
    foo = foo)
}
