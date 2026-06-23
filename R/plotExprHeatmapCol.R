# plotExprHeatmapCol.R  — CATALYST-free rewrite of plotExprHeatmapCol.
# Public name: plotExprHeatmapCol (at cut-over, a deprecated `plotExprHeatmapCol`
# alias will forward here with a .Deprecated() warning for backward compatibility).
# -----------------------------------------------------------------------------
# The function BODY is David Priest's own code (extracted from CustomFunctions);
# it is not derived from CATALYST. Only its calls into CATALYST internals are
# replaced here with the MIT .wl_* equivalents, and the
#   environment(...) <- asNamespace('CATALYST')
# hack is removed. External generics are namespace-qualified so the function
# resolves with no CATALYST on the search path.
#
# Requires (already seekit deps): ComplexHeatmap, circlize, grid,
# RColorBrewer, magrittr, SummarizedExperiment, and the winglab helpers
# subSCE / anno_clusters2 / anno_counts1 (all CATALYST-free) + internals.R.
# -----------------------------------------------------------------------------

plotExprHeatmapCol <- function (x, features = NULL, by = c("sample_id", "cluster_id",
    "both"), k = "meta20", assay = "exprs", fun = c("median",
    "mean", "sum"), scale = c("first", "last", "never"), q = 0.01, sub = FALSE,
    row_anno = TRUE, row_names = FALSE, col_anno = TRUE, row_clust = TRUE, col_clust = TRUE,
    row_dend = TRUE, col_dend = TRUE, bars = FALSE, perc = FALSE, title = " ",
    bin_anno = FALSE, hm_pal = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
    k_pal = NULL, named_k_pal = NULL, distance = c("euclidean",
        "maximum", "manhattan", "canberra", "binary", "minkowski"),
    linkage = c("average", "ward.D", "single", "complete", "mcquitty",
        "median", "centroid", "ward.D2"), row_order = NULL)
{
    if (sub == TRUE)
    {
      x <- subSCE(x, n_cells = 1000)
      title <- "WARNING: downsampled SCE used."
    }

    args <- as.list(environment())
    distance <- match.arg(distance)
    linkage <- match.arg(linkage)
    scale <- match.arg(scale)
    fun <- match.arg(fun)
    by <- match.arg(by)
    x <- x[unique(.wl_get_features(x, features)), ]

    if (by != "sample_id") {
        x$cluster_id <- .wl_resolve_cluster_ids(x, k)   # [#7] shared resolver
    }
    if (by == "both")
        by <- c("cluster_id", "sample_id")

    .do_agg <- function() {
        z <- .wl_agg(x, by, fun, assay)
        if (length(by) == 1)
            return(z)
        magrittr::set_rownames(do.call("rbind", z), levels(x$cluster_id))
    }
    .do_scale <- function() {
        if (scale == "first") {
            z <- SummarizedExperiment::assay(x, assay)
            z <- .wl_scale_exprs(z, 1, q)
            SummarizedExperiment::assay(x, assay, withDimnames = FALSE) <- z
            return(x)
        }
        else .wl_scale_exprs(z, 1, q)
    }
    z <- switch(scale, first = {
        x <- .do_scale()
        .do_agg()
    }, last = {
        z <- .do_agg()
        .do_scale()
    }, never = {
        .do_agg()
    })

    if (length(by) == 1)
        z <- t(z)
    if (scale != "never" && !(assay == "counts" && fun == "sum")) {
        qs <- round(stats::quantile(z, c(0.01, 0.99)) * 5)/5
        lgd_aes <- list(at = seq(qs[1], qs[2], 0.2))
    }
    else lgd_aes <- list()
    lgd_aes$title_gp <- grid::gpar(fontsize = 10, fontface = "bold",
        lineheight = 0.8)

    # Apply custom row order
    if (!is.null(row_order)) {
      z <- z[row_order, , drop = FALSE]
    }

    # Extract the cluster IDs in the order of the heatmap rows
    cluster_order <- rownames(z)

    # Ensure named_k_pal is provided and matches the cluster_order
    if (!is.null(named_k_pal)) {
        k_pal <- named_k_pal[cluster_order]
    } else {
        stop("Please provide a named color palette 'named_k_pal' matching your clusters.")
    }

    # Create the left annotation
    if (!isFALSE(row_anno)) {
        left_anno <- anno_clusters2(cluster_order, k_pal)
    } else {
        left_anno <- NULL
    }

    if (!isFALSE(col_anno) && length(by) == 2) {
        top_anno <- .wl_anno_factors(x, levels(x$sample_id), col_anno,
            "column")
    } else {
        top_anno <- NULL
    }

    if (bars) {
        right_anno <- anno_counts1(x[[by[1]]], perc)
    } else {
        right_anno <- NULL
    }

    if (bin_anno) {
        cell_fun <- function(j, i, x, y, ...) grid::grid.text(gp = grid::gpar(fontsize = 6),
            sprintf("%.2f", z[i, j]), x, y)
    } else {
        cell_fun <- NULL
    }

    a <- ifelse(assay == "exprs", "expression", assay)
    f <- switch(fun, median = "med", fun)
    hm_title <- switch(scale, first = sprintf("%s %s\n%s", fun,
        "scaled", a), last = sprintf("%s %s\n%s", "scaled", fun,
        a), never = paste(fun, a, sep = "\n"))
    if (length(by) == 2) {
        col_title <- features
    } else if (length(features) == 1 && features %in% c("type",
        "state")) {
        col_title <- paste0(features, "_markers")
    } else {
        col_title <- ""
    }

    p <- ComplexHeatmap::Heatmap(matrix = z, name = hm_title,
        col = circlize::colorRamp2(seq(min(z),
        max(z), length.out = 100),
        grDevices::colorRampPalette(hm_pal)(100)),
        column_title = col_title, column_title_side = ifelse(length(by) == 2, "top", "bottom"),
        cell_fun = cell_fun,
        cluster_rows = row_clust,
        cluster_columns = col_clust,
        show_row_dend = row_dend,
        show_column_dend = col_dend,
        clustering_distance_rows = distance,
        clustering_method_rows = linkage,
        clustering_distance_columns = distance,
        clustering_method_columns = linkage,
        show_row_names = row_names,
        row_names_side = ifelse(by[1] == "cluster_id" || isFALSE(row_anno) && !row_dend || isFALSE(row_clust), "left", "right"),
        top_annotation = top_anno,
        left_annotation = left_anno,
        right_annotation = right_anno,
        rect_gp = grid::gpar(col = "white"),
        heatmap_legend_param = lgd_aes,
        row_names_max_width = grid::unit(13, "cm"),
        row_title = title,
        column_names_rot = 45,
        column_names_gp = grid::gpar(fontsize = 12),
        row_names_gp = grid::gpar(fontsize = 14))

    ht <- ComplexHeatmap::draw(p)

    row_labels_after_clustering <- rownames(z)[ComplexHeatmap::row_order(ht)]
    row_labels_df <- data.frame(row_labels = row_labels_after_clustering)

    list(heatmap = p, row_labels_df = row_labels_df)
}
