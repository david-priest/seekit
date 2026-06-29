# plotMultiExprHeatmap.R — CATALYST-free rewrite of the CustomFunctionsSept25
# function (my_functions.R). Two side-by-side expression heatmaps sharing one
# colour scale and one row (cluster) order: hm1 on the left (with the cluster
# colour annotation), hm2 on the right (markers, optionally with a count bar).
#
# Clean-room changes vs the legacy version:
#   - CATALYST:::.get_features / .agg / .scale_exprs  ->  .wl_get_features /
#     .wl_agg / .wl_scale_exprs (R/wl_internals.R). Same documented behaviour.
#   - anno_clusters2 / anno_counts1 are the package's CATALYST-free annotations
#     (R/heatmap_annotations.R).
#   - The `environment(...) <- asNamespace('CATALYST')` hack is removed.
#   - hm2 = "abundances" used CATALYST::plotFreqHeatmap and is NOT supported here
#     (a clear error is raised). Pass marker names for `hm2` instead — which is
#     what every CMV/AIM call already does. 2026-06 seekit migration.

plotMultiExprHeatmap <- function(x,
                                 hm1 = "type",
                                 hm2 = "abundances",
                                 k = "meta20",
                                 m = NULL,
                                 assay = "exprs",
                                 fun = "median",
                                 scale = "first",
                                 q = 0.01,
                                 named_k_pal = NULL,
                                 row_clust = TRUE,
                                 col_clust = c(TRUE, TRUE),
                                 row_dend = TRUE,
                                 col_dend = TRUE,
                                 row_names = FALSE,
                                 row_names_side = "left",
                                 show_column_titles = TRUE,
                                 row_order = NULL,
                                 bars = TRUE,
                                 perc = TRUE,
                                 hm_pal = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                                 ...) {

    # --- Argument Validation ---
    if (is.null(named_k_pal)) {
        stop("Please provide a named color palette 'named_k_pal' matching your clusters.")
    }
    if (length(col_clust) == 1) col_clust <- rep(col_clust, 2)
    if (length(col_dend) == 1) col_dend <- rep(col_dend, 2)

    # Internal helper function to aggregate data
    .get_agg_matrix <- function(features_in, by = "cluster_id") {
        sce_sub <- x[unique(.wl_get_features(x, features_in)), ]
        if (by != "sample_id") sce_sub$cluster_id <- sce_sub[[k]]

        z <- .wl_agg(sce_sub, by, fun, assay)

        if (scale == "first") {
            z_scaled <- assay(sce_sub, assay)
            z_scaled <- .wl_scale_exprs(z_scaled, 1, q)
            assay(sce_sub, assay, FALSE) <- z_scaled
            z <- .wl_agg(sce_sub, by, fun, assay)
        } else if (scale == "last") {
            z <- .wl_scale_exprs(z, 1, q)
        }

        return(t(z))
    }

    # --- Data Preparation ---
    z1 <- NULL
    z2 <- NULL

    if (!isFALSE(hm1)) {
        z1 <- .get_agg_matrix(hm1)
    }

    is_hm2_expr <- !(isTRUE(hm2 == "abundances"))
    if (is_hm2_expr) {
        features_hm2 <- if(isTRUE(hm2 == "state")) "state" else hm2
        z2 <- .get_agg_matrix(features_hm2)
    }

    # --- Unified Color Scale Calculation ---
    color_fun <- NULL
    if (!is.null(z1) && !is.null(z2)) {
        global_min <- min(c(as.matrix(z1), as.matrix(z2)), na.rm = TRUE)
        global_max <- max(c(as.matrix(z1), as.matrix(z2)), na.rm = TRUE)
        message(paste("Unified scale created with range:", round(global_min, 2), "to", round(global_max, 2)))
        color_fun <- colorRamp2(seq(global_min, global_max, length.out = 100),
                                colorRampPalette(hm_pal)(100))
    }

    # --- Heatmap Generation ---
    ht_list <- NULL

    # --- Generate Heatmap 1 ---
    if (!is.null(z1)) {
        if (row_clust && is.null(row_order)) {
            temp_ht <- Heatmap(z1, cluster_rows = TRUE)
            ht_drawn <- draw(temp_ht, merge_legend = TRUE, newpage = FALSE)
            row_order <- row_order(ht_drawn)
            message("Capturing row order from the first heatmap to sync subsequent heatmaps.")
        }

        if (!is.null(row_order)) {
            z1 <- z1[row_order, , drop = FALSE]
            row_clust <- FALSE
        }

        cluster_order <- rownames(z1)
        k_pal_ordered <- named_k_pal[cluster_order]
        left_anno <- anno_clusters2(cluster_order, k_pal_ordered)

        a <- ifelse(assay == "exprs", "expression", assay)
        hm_title <- switch(scale,
                           first = sprintf("%s scaled %s", fun, a),
                           last = sprintf("scaled %s %s", fun, a),
                           never = paste(fun, a, sep = "\n"))

        ht1 <- Heatmap(
            matrix = z1,
            name = hm_title,
            col = color_fun,
            cluster_rows = row_clust,
            show_row_dend = row_dend && row_clust,
            cluster_columns = col_clust[1],
            show_column_dend = col_dend[1],
            show_row_names = row_names,
            row_names_side = row_names_side,
            left_annotation = left_anno,
            column_title = if (show_column_titles) paste(hm1, collapse = ", ") else NULL,
            rect_gp = gpar(col = "white"),
            heatmap_legend_param = list(title = hm_title),
            show_heatmap_legend = TRUE
        )
        ht_list <- ht1
    }

    # --- Generate Heatmap 2 ---
    if (is_hm2_expr && !is.null(z2)) {
        if (!is.null(row_order)) {
            z2 <- z2[row_order, , drop = FALSE]
        }

        right_anno <- if (bars) anno_counts1(x[[k]], perc) else NULL

        ht2 <- Heatmap(
            matrix = z2,
            name = " ",
            col = color_fun,
            cluster_rows = FALSE,
            show_row_dend = FALSE,
            cluster_columns = col_clust[2],
            show_column_dend = col_dend[2],
            show_row_names = FALSE,
            right_annotation = right_anno,
            column_title = if (show_column_titles) paste(hm2, collapse = ", ") else NULL,
            rect_gp = gpar(col = "white"),
            show_heatmap_legend = FALSE
        )
        ht_list <- ht_list + ht2

    } else if (!is_hm2_expr) {
        stop("plotMultiExprHeatmap(): hm2 = \"abundances\" is not supported in the ",
             "CATALYST-free seekit build (it relied on CATALYST::plotFreqHeatmap). ",
             "Pass marker names for `hm2` instead (e.g. hm2 = c(\"Ki67\", \"CD25\")), ",
             "or use the abundance plotters (plotAbundanceStacked / plotAbundancesDiff).",
             call. = FALSE)
    }

    return(ht_list)
}
