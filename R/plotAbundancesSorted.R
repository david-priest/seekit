# plotAbundancesSorted — migrated from CyTOF nBass_helpers.R into seekit (CATALYST-free).
# 2026-06-10: lifted verbatim, de-CATALYST'd (.wl_* internals, namespace hack removed).

plotAbundancesSorted <- function (x, k = "meta20", by = c("sample_id", "cluster_id"),
                                  group_by = "condition", shape_by = NULL,
                                  shape_palette = NULL,
                                  point_fill = "red", point_color = "black",
                                  col_clust = TRUE, n_cols = 4, log = FALSE,
                                  miny = 0.01, maxy = NA, point_size = 3,
                                  distance = c("euclidean", "maximum", "manhattan",
                                               "canberra", "binary", "minkowski"),
                                  linkage = c("average", "ward.D", "single", "complete",
                                              "mcquitty", "median", "centroid", "ward.D2"),
                                  k_pal = .wl_cluster_cols, clusters_order = NULL) {
  by <- match.arg(by)
  .wl_check_sce(x, TRUE)
  k <- .wl_check_k(x, k)
  .wl_check_cd_factor(x, group_by)
  .wl_check_cd_factor(x, shape_by)
  .wl_check_pal(k_pal)
  linkage <- match.arg(linkage)
  distance <- match.arg(distance)
  stopifnot(is.logical(col_clust), length(col_clust) == 1)

  shapes <- .wl_get_shapes(x, shape_by)
  if (is.null(shapes)) shape_by <- NULL

  if (by == "sample_id") {
    nk <- nlevels(.wl_cluster_ids(x, k))
    if (length(k_pal) < nk) k_pal <- colorRampPalette(k_pal)(nk)
  }

  ns <- table(cluster_id = .wl_cluster_ids(x, k), sample_id = .wl_sample_ids(x))
  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)
  m <- match(df$sample_id, x$sample_id)
  for (i in c(shape_by, group_by)) df[[i]] <- x[[i]][m]

  if (by == "sample_id" && col_clust && length(unique(df$sample_id)) > 1) {
    d <- dist(t(fq), distance); h <- hclust(d, linkage)
    df$sample_id <- factor(df$sample_id, colnames(fq)[h$order])
  }
  if (!is.null(clusters_order)) df$cluster_id <- factor(df$cluster_id, levels = clusters_order)
  if (log == TRUE) df$Freq <- df$Freq + 0.02

  p <- ggplot(df, aes_string(y = "Freq")) +
    labs(x = NULL, y = "Proportion [%]") + theme_bw() +
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"),
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = NA, color = NA),
          axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
          legend.key.height = unit(0.8, "lines"))

  p <- p + ggh4x::facet_wrap2(~cluster_id, scales = "free_y", ncol = n_cols, axes = "all") +
    geom_boxplot(aes_string(x = group_by), color = "grey16", fill = "grey90",
                 position = position_dodge(), alpha = 0.8,
                 outlier.color = NA, show.legend = FALSE)

  if (!is.null(shape_by)) {
    p <- p + geom_quasirandom(aes_string(x = group_by, shape = shape_by),
                              size = point_size, width = 0.2,
                              fill = point_fill, color = point_color)
    if (!is.null(shape_palette)) p <- p + scale_shape_manual(values = shape_palette)
  } else {
    p <- p + geom_quasirandom(aes_string(x = group_by), fill = "grey84",
                              size = point_size, width = 0.2, shape = 21)
  }

  if (log == TRUE) {
    p + scale_y_continuous(trans = "log10", limits = c(miny, maxy),
                           breaks = c(0.01, 0.1, 1, 10, 100), labels = c(0.01, 0.1, 1, 10, 100)) +
      annotation_logticks(base = 10, sides = "l", outside = TRUE) +
      coord_cartesian(clip = "off") + scale_size_area(max_size = 15) +
      theme(axis.text.y = element_text(margin = margin(r = 8)))
  } else {
    p + scale_y_continuous(limits = c(0, maxy)) +
      coord_cartesian(clip = "off") + scale_size_area(max_size = 15)
  }
}
