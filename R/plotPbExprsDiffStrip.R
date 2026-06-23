# plotPbExprsDiffStrip
#
# Sibling of plotPbExprsDiff with a TRUE separated stats strip per cell:
# each (cluster, antigen) facet is composed of two stacked sub-plots,
#
#   [ stats strip ]   <- brackets only, no axes, fixed-height across the grid
#   [   data plot ]   <- boxplots + points, axes intact, untouched by brackets
#
# Stacked via patchwork::plot_layout(heights = ...) per cell, then the cells
# are arranged in a 2D grid (rows = cluster_id, cols = antigen by default;
# swap = TRUE flips them). Modeled on plotAbundancesDiffStrip.R.
#
# Grid-label handling (because patchwork cells don't share strip text like
# facet_grid2 does):
#   - Top-row cells get the column variable name as the strip title.
#   - Leftmost-column cells get the row variable name in the y-axis title.
#   - Interior cells suppress both labels, so the visual layout matches
#     facet_grid2 with strip.placement = "outside".
#
# Properties (mirrors plotAbundancesDiffStrip):
#   - Strip y-axis is 0 .. max_brackets across ALL cells -> identical bracket
#     spacing throughout the grid; cells with fewer brackets place them at
#     the bottom of the strip (nearest the data plot below).
#   - Data plot panel width is locked in cm; aspect ratio sets the height.
#   - tight_top_axis terminates each cell's y-axis at the next pretty tick
#     above its max data point.
#   - hide_ns is a *render* knob; the full unfiltered stats are stashed on
#     the returned patchwork object as attr(p, "source_stats") so f2() can
#     export exact p-values for every comparison (including ns).
#
# Dependencies (loaded by SETUP 1 of the calling .Rmd):
#   SingleCellExperiment, ggplot2, ggbeeswarm, patchwork, ggh4x, scales, grid,
#   rstatix, dplyr, rlang, reshape2 (for melt, via CATALYST namespace).
# Optional (only when use_richtext = TRUE): ggtext.

plotPbExprsDiffStrip <- function(
    x,
    k                  = "meta20",
    features           = "state",
    assay              = "exprs",
    mean_or_med        = c("median", "mean", "sum"),
    color_by           = "condition",
    group_by           = color_by,            # for the x-axis / stats grouping
    group_levels       = NULL,
    clusters_order     = NULL,
    features_order     = NULL,                # explicit antigen order; defaults to `features`
    excluded_clusters  = NULL,
    # Whitelist filter applied AFTER the pseudo-bulk per-(cluster, sample)
    # expression values are computed, so the values themselves are unaffected
    # -- only the SET of clusters whose rows get rendered is reduced. External
    # stats passed via `external_stats` are auto-restricted to the kept
    # clusters by the existing per-cell stats-prep block (no change needed
    # downstream). Use case: paper-figure subsets (e.g. "only show CD4 CTL"
    # or "only TEMRA + EM Ki67+"). NULL keeps all clusters (current
    # behaviour). If supplied AND `clusters_order` is NULL, the order of
    # `keep_clusters` is used as the row order; pass `clusters_order`
    # explicitly to override.
    keep_clusters      = NULL,
    fill_palette       = NULL,
    point_size         = 2,
    facet_ratio        = 1,
    panel_width_cm     = 3,
    swap               = FALSE,               # FALSE = cluster_id rows, antigen cols
    merging_col        = NULL,                # mirrors plotPbExprsDiff; resolves k from colData
    # Stats
    external_stats     = NULL,
    show_stats         = TRUE,
    hide_ns            = TRUE,
    # Multiple-testing correction view. Three options:
    #   "per_pair"    -> diffcyt default; BH across clusters within each pair
    #                   (m = clusters in that pair-call). Uses p.adj.signif.
    #   "per_cluster" -> recommended (per PI); BH across pairs within each
    #                   (cluster, marker) only. Doesn't correct across clusters
    #                   or markers. Uses p.adj_per_cluster.signif.
    #   "global"      -> BH across the full clusters × markers × pairs family.
    #                   Uses p.adj_global.signif.
    # Convenience switch over `label_col`; if `label_col` is set explicitly
    # to anything other than the default, it wins.
    correction         = c("per_pair", "per_cluster", "global"),
    # Asterisk / label styling
    label_col          = "p.adj.signif",
    label_size         = 3,
    label_fontface     = "plain",
    use_richtext       = TRUE,
    asterisk_pt_multiplier  = 1.4,
    asterisk_vjust     = 0.3,
    label_vjust        = 0,
    asterisk_y_offset  = 0,
    bracket_color      = "black",
    bracket_size       = 0.4,
    bracket_tip_npc    = 0.20,
    label_nudge_npc    = 0.05,
    # Strip layout
    strip_height_fraction = 0.25,
    panel_gap_pt          = 8,
    # Y-axis
    y_expand_low_mult  = 0.02,
    y_expand_high_mult = 0.05,
    tight_top_axis           = TRUE,
    tight_top_axis_overhang  = 0.02,
    # Title placement / typography
    title_above_strip  = TRUE,
    title_size         = 10,
    axis_text_size     = 8,
    row_label_size     = NULL,                # defaults to title_size if NULL
    # Theme
    nature_style       = TRUE,
    legend_position    = "right",
    # y-axis label root (the original showed "median expression" / "mean expression")
    y_axis_root        = NULL,
    # Diagnostic
    verbose            = TRUE
) {
    # ---- Package dependencies --------------------------------------------
    if (!requireNamespace("patchwork", quietly = TRUE))
        stop("patchwork is required: install.packages('patchwork').")
    if (!requireNamespace("ggh4x", quietly = TRUE))
        stop("ggh4x is required (force_panelsizes): install.packages('ggh4x').")
    if (!requireNamespace("reshape2", quietly = TRUE))
        stop("reshape2 is required (melt): install.packages('reshape2').")
    if (isTRUE(use_richtext) && !requireNamespace("ggtext", quietly = TRUE))
        stop("ggtext is required when use_richtext = TRUE: install.packages('ggtext').")

    mean_or_med <- match.arg(mean_or_med)
    if (is.null(row_label_size)) row_label_size <- title_size

    # Correction-view convenience switch (see param doc above).
    correction <- match.arg(correction)
    if (identical(label_col, "p.adj.signif")) {
        label_col <- switch(correction,
                            per_pair    = "p.adj.signif",
                            per_cluster = "p.adj_per_cluster.signif",
                            global      = "p.adj_global.signif")
    }

    # ---- Cluster ids (mirrors plotPbExprsDiff logic) ----------------------
    if (!is.null(merging_col) || k %in% names(colData(x))) {
        cl_vec <- factor(x[[k]])
    } else {
        .wl_check_sce(x)
        k <- .wl_check_k(x, k)
        cl_vec <- .wl_cluster_ids(x, k)
    }

    # ---- Sanity checks ----------------------------------------------------
    .wl_check_assay(x, assay)
    .wl_check_cd_factor(x, color_by)
    .wl_check_cd_factor(x, group_by)

    # ---- Build long-format df of per-sample-cluster-antigen values --------
    x <- x[.wl_get_features(x, features), ]
    x$cluster_id <- cl_vec
    by <- c("cluster_id", "sample_id")
    ms <- .wl_agg(x, by, mean_or_med, assay)
    df <- reshape2::melt(ms, varnames = c("antigen", by[length(by)]))

    if (!is.null(cl_vec)) {
        df$cluster_id <- df$L1
    }

    i <- match(df$sample_id, x$sample_id)
    j <- setdiff(names(colData(x)), c(names(df), "cluster_id"))
    df <- cbind(df, colData(x)[i, j, drop = FALSE])

    ncs <- table(as.list(colData(x)[by]))
    ncs <- rep(c(t(ncs)), each = nrow(x))
    df  <- df[ncs > 0, , drop = FALSE]

    if (!is.null(excluded_clusters)) {
        df <- df[!df$cluster_id %in% excluded_clusters, ]
    }
    if (!is.null(keep_clusters)) {
        # Whitelist filter. Applied after the pseudo-bulk values are computed,
        # so per-cluster expression values are unaffected -- only the rendered
        # subset changes. external_stats auto-narrows downstream because the
        # stats-prep block keys on `clusters` (the post-filter levels).
        .missing <- setdiff(keep_clusters, unique(df$cluster_id))
        if (length(.missing) > 0) {
            warning("keep_clusters contains cluster ids not present in data: ",
                    paste(.missing, collapse = ", "))
        }
        df <- df[df$cluster_id %in% keep_clusters, ]
    }

    # ---- Order clusters + antigens ----------------------------------------
    if (!is.null(clusters_order)) {
        df$cluster_id <- factor(df$cluster_id, levels = clusters_order)
    } else if (!is.null(keep_clusters)) {
        # No explicit order but keep_clusters was supplied -- use its order
        # as the row layout (intersect drops names not present in data).
        df$cluster_id <- factor(
            df$cluster_id,
            levels = intersect(keep_clusters, unique(as.character(df$cluster_id)))
        )
    } else {
        df$cluster_id <- factor(df$cluster_id, levels = unique(as.character(df$cluster_id)))
    }
    df$cluster_id <- droplevels(df$cluster_id)
    clusters <- levels(df$cluster_id)

    if (!is.null(features_order)) {
        df$antigen <- factor(df$antigen, levels = features_order)
    } else {
        # `features` may be a marker class shortcut ("state" / "type"); use whatever
        # ended up in the long-format df.
        df$antigen <- factor(df$antigen, levels = unique(as.character(df$antigen)))
    }
    df$antigen <- droplevels(df$antigen)
    antigens <- levels(df$antigen)

    # ---- Group levels on x-axis -------------------------------------------
    if (is.null(group_levels)) {
        group_levels <- if (is.factor(df[[group_by]])) {
            levels(droplevels(factor(df[[group_by]])))
        } else {
            sort(unique(as.character(df[[group_by]])))
        }
    }
    df[[group_by]] <- factor(df[[group_by]], levels = group_levels)

    # ---- n samples per condition (Immunity figure legend reporting) -------
    # Counts distinct sample_ids per group_by level after all filters.
    n_per_cond_df <- df %>%
        dplyr::distinct(sample_id, !!rlang::sym(group_by)) %>%
        dplyr::count(!!rlang::sym(group_by), name = "n_samples") %>%
        as.data.frame()
    message("plotPbExprsDiffStrip: n samples per ", group_by, ":")
    print(n_per_cond_df, row.names = FALSE)

    # ---- Y-axis root label -----------------------------------------------
    if (is.null(y_axis_root)) {
        y_axis_root <- paste(mean_or_med, if (assay == "exprs") "expression" else assay)
    }

    # ---- Pre-scan per-cell stats + global max_brackets -------------------
    # cell_stats[[cluster, antigen]] holds the rendered (post-hide_ns) rows.
    cell_stats <- list()
    full_stats_rows <- list()
    if (isTRUE(show_stats) && !is.null(external_stats)) {
        es <- as.data.frame(external_stats)
        # Restrict the full stash to clusters/antigens actually plotted.
        keep <- es$cluster_id %in% clusters & es$antigen %in% antigens
        es_plot <- es[keep, , drop = FALSE]
        full_stats_rows[[1]] <- es_plot
        for (cl in clusters) for (ag in antigens) {
            sub <- es_plot[as.character(es_plot$cluster_id) == cl &
                           as.character(es_plot$antigen)    == ag, , drop = FALSE]
            if (isTRUE(hide_ns)) {
                sub <- sub[as.character(sub[[label_col]]) != "ns", , drop = FALSE]
            }
            cell_stats[[paste(cl, ag, sep = "\037")]] <- sub
        }
    }
    max_brackets <- 0L
    if (length(cell_stats)) {
        max_brackets <- max(0L, vapply(cell_stats,
                                       function(s) if (is.null(s)) 0L else nrow(s),
                                       integer(1)))
    }

    # ---- Strip / data height fractions (uniform across cells) ------------
    strip_h <- strip_height_fraction
    data_h  <- 1 - strip_height_fraction
    strip_aspect <- facet_ratio * strip_h / data_h

    # ---- Label transform (asterisk styling, copied from plotAbundancesDiffStrip) ----
    if (isTRUE(use_richtext)) {
        asterisk_pt <- max(1, round(label_size * 3.5 * asterisk_pt_multiplier))
        label_transform <- function(s) {
            has_ast <- grepl("[*]", s)
            if (any(has_ast)) {
                s_esc <- gsub("\\*", "&#42;", s)
                s[has_ast] <- paste0(
                    "<span style='font-size:", asterisk_pt, "pt; line-height:0.6'>",
                    s_esc[has_ast], "</span>"
                )
            }
            s
        }
    } else {
        label_transform <- function(s) s
    }

    # ---- Per-cell builders ------------------------------------------------
    # show_row_label / show_col_title control which cells own the row label
    # (cluster on left) and column title (antigen on top). In swap = FALSE
    # mode rows = clusters, cols = antigens; in swap = TRUE rows = antigens,
    # cols = clusters.
    build_data_plot <- function(df_cell, row_lbl, col_lbl,
                                show_row_label, show_col_title) {
        y_lab <- if (isTRUE(show_row_label)) paste0(row_lbl, "\n", y_axis_root) else NULL
        p <- ggplot2::ggplot(
                df_cell,
                ggplot2::aes(x    = .data[[group_by]],
                             y    = value,
                             fill = .data[[color_by]])
            ) +
            ggplot2::geom_boxplot(
                color = "black", linewidth = 0.3, width = 0.75,
                alpha = 0.9, outlier.color = NA,
                show.legend = TRUE
            ) +
            ggbeeswarm::geom_quasirandom(
                shape = 21, fill = "grey84", color = "black", stroke = 0.5,
                size = point_size, width = 0.2, alpha = 0.8
            ) +
            ggplot2::labs(
                x = NULL,
                y = y_lab,
                title = if (isTRUE(show_col_title) && !isTRUE(title_above_strip))
                            col_lbl else NULL
            ) +
            (if (isTRUE(tight_top_axis)) {
                y_max <- max(df_cell$value, na.rm = TRUE)
                if (!is.finite(y_max) || y_max <= 0) y_max <- 1
                breaks_vec <- pretty(c(0, y_max * 1.001), n = 5)
                top_break  <- max(breaks_vec)
                axis_top   <- top_break * (1 + tight_top_axis_overhang)
                ggplot2::scale_y_continuous(
                    breaks = breaks_vec,
                    limits = c(0, axis_top),
                    expand = ggplot2::expansion(mult = c(y_expand_low_mult, 0))
                )
            } else {
                ggplot2::scale_y_continuous(
                    limits = c(0, NA),
                    expand = ggplot2::expansion(mult = c(y_expand_low_mult, y_expand_high_mult))
                )
            }) +
            ggplot2::scale_x_discrete(limits = group_levels, drop = FALSE) +
            ggh4x::force_panelsizes(cols = grid::unit(panel_width_cm, "cm")) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                panel.grid       = ggplot2::element_blank(),
                aspect.ratio     = facet_ratio,
                plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5,
                                                         size = title_size,
                                                         margin = ggplot2::margin(b = 2)),
                axis.title.y     = ggplot2::element_text(size = row_label_size,
                                                         face = "bold"),
                axis.text        = ggplot2::element_text(color = "black", size = axis_text_size),
                axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1,
                                                         size = axis_text_size),
                axis.ticks       = ggplot2::element_line(color = "black"),
                plot.margin      = ggplot2::margin(2, 2, 2, 2)
            )

        if (!is.null(fill_palette)) {
            p <- p + ggplot2::scale_fill_manual(values = fill_palette, name = color_by, drop = FALSE)
        }

        if (isTRUE(nature_style)) {
            p <- p + ggplot2::theme(
                panel.border     = ggplot2::element_blank(),
                panel.background = ggplot2::element_blank(),
                axis.line        = ggplot2::element_line(color = "black", linewidth = 0.4)
            )
        } else {
            p <- p + ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
            )
        }
        p
    }

    build_strip_plot <- function(stat_cell, col_lbl, show_col_title) {
        n_brackets <- if (is.null(stat_cell)) 0 else nrow(stat_cell)
        n_groups   <- length(group_levels)
        y_axis_max <- max(max_brackets, 1)

        if (n_brackets > 0) {
            bracket_y <- n_brackets - seq_len(n_brackets) + 0.5
            x1 <- match(as.character(stat_cell$group1), group_levels)
            x2 <- match(as.character(stat_cell$group2), group_levels)
            x_mid <- (x1 + x2) / 2
            raw_labels <- as.character(stat_cell[[label_col]])
            labels     <- label_transform(raw_labels)

            tip_dy   <- bracket_tip_npc
            nudge_dy <- label_nudge_npc

            seg_h <- data.frame(x = x1, xend = x2, y = bracket_y, yend = bracket_y)
            seg_l <- data.frame(x = x1, xend = x1, y = bracket_y, yend = bracket_y - tip_dy)
            seg_r <- data.frame(x = x2, xend = x2, y = bracket_y, yend = bracket_y - tip_dy)
            seg_df <- rbind(seg_h, seg_l, seg_r)

            is_ast <- grepl("[*]", raw_labels)
            label_y <- bracket_y + nudge_dy
            label_y[is_ast] <- label_y[is_ast] + asterisk_y_offset
            vj_vec <- ifelse(is_ast, asterisk_vjust, label_vjust)
            text_df <- data.frame(x = x_mid, y = label_y, label = labels, vj = vj_vec)
        } else {
            seg_df  <- data.frame(x = numeric(0), xend = numeric(0),
                                  y = numeric(0), yend = numeric(0))
            text_df <- data.frame(x = numeric(0), y = numeric(0),
                                  label = character(0), vj = numeric(0))
        }

        p <- ggplot2::ggplot()
        if (nrow(seg_df) > 0) {
            # lineend = "round" rounds the bottoms of the vertical tips AND
            # makes the horizontal bar's caps meet the verticals without a
            # corner gap. Same convention as plotAbundancesDiffStrip.
            p <- p + ggplot2::geom_segment(
                data    = seg_df,
                mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                color   = bracket_color, linewidth = bracket_size,
                lineend = "round", linejoin = "round"
            )
            if (isTRUE(use_richtext)) {
                p <- p + ggtext::geom_richtext(
                    data    = text_df,
                    mapping = ggplot2::aes(x = x, y = y, label = label, vjust = vj),
                    size    = label_size, fontface = label_fontface,
                    label.size    = NA, fill = NA,
                    label.padding = grid::unit(c(0, 0, 0, 0), "lines")
                )
            } else {
                p <- p + ggplot2::geom_text(
                    data    = text_df,
                    mapping = ggplot2::aes(x = x, y = y, label = label, vjust = vj),
                    size    = label_size, fontface = label_fontface
                )
            }
        }

        p <- p +
            ggplot2::labs(title = if (isTRUE(show_col_title) && isTRUE(title_above_strip))
                                       col_lbl else NULL) +
            ggplot2::scale_x_continuous(limits = c(0.5, n_groups + 0.5), expand = c(0, 0)) +
            ggplot2::scale_y_continuous(limits = c(0, y_axis_max + 0.5),
                                        expand = c(0, 0),
                                        oob    = scales::oob_keep) +
            ggplot2::coord_cartesian(clip = "off") +
            ggh4x::force_panelsizes(cols = grid::unit(panel_width_cm, "cm")) +
            ggplot2::theme_void() +
            ggplot2::theme(
                aspect.ratio = strip_aspect,
                plot.title   = ggplot2::element_text(face = "bold", hjust = 0.5,
                                                     size = title_size,
                                                     margin = ggplot2::margin(b = 2)),
                plot.margin  = ggplot2::margin(0, 2, 0, 2)
            )
        p
    }

    # ---- Build grid of composites ----------------------------------------
    # In swap = FALSE: rows = clusters, cols = antigens. ncol = length(antigens).
    # In swap = TRUE:  rows = antigens, cols = clusters. ncol = length(clusters).
    if (!isTRUE(swap)) {
        row_vals <- clusters; col_vals <- antigens
        get_row <- function(cl, ag) cl
        get_col <- function(cl, ag) ag
    } else {
        row_vals <- antigens; col_vals <- clusters
        get_row <- function(cl, ag) ag
        get_col <- function(cl, ag) cl
    }
    n_cols_grid <- length(col_vals)

    composites <- list()
    idx <- 0L
    for (rv in row_vals) for (cv in col_vals) {
        idx <- idx + 1L
        # Resolve which is cluster vs antigen for data extraction
        if (!isTRUE(swap)) { cl <- rv; ag <- cv } else { cl <- cv; ag <- rv }

        df_cell <- df[as.character(df$cluster_id) == cl &
                      as.character(df$antigen)    == ag, , drop = FALSE]
        stat_cell <- cell_stats[[paste(cl, ag, sep = "\037")]]

        # Position in grid (1-indexed)
        row_pos <- ((idx - 1L) %/% n_cols_grid) + 1L
        col_pos <- ((idx - 1L) %%  n_cols_grid) + 1L
        show_col_title  <- (row_pos == 1L)
        show_row_label  <- (col_pos == 1L)

        # Skip cells with no data entirely? Keep them for layout consistency
        # but they'll just render an empty axis. Set y-range manually so the
        # plot doesn't error.
        if (nrow(df_cell) == 0) {
            df_cell <- data.frame(value = c(0, 1))
            df_cell[[group_by]] <- factor(group_levels[c(1, 1)], levels = group_levels)
            df_cell[[color_by]] <- factor(group_levels[c(1, 1)], levels = group_levels)
        }

        p_data  <- build_data_plot(df_cell, row_lbl = rv, col_lbl = cv,
                                   show_row_label = show_row_label,
                                   show_col_title = show_col_title)
        p_strip <- build_strip_plot(stat_cell, col_lbl = cv,
                                    show_col_title = show_col_title)

        composites[[idx]] <- p_strip / p_data +
            patchwork::plot_layout(heights = c(strip_h, data_h))
    }

    # ---- Combine via wrap_plots ------------------------------------------
    gap <- panel_gap_pt
    out <- patchwork::wrap_plots(composites, ncol = n_cols_grid) &
        ggplot2::theme(
            legend.position = legend_position,
            plot.margin     = ggplot2::margin(t = gap, r = gap, b = gap, l = gap, unit = "pt")
        )
    if (legend_position != "none") {
        out <- out + patchwork::plot_layout(guides = "collect")
    }

    # ---- Stash source data + stats for f2() Excel export ----------------
    attr(out, "source_data") <- df
    attr(out, "source_stats") <- if (length(full_stats_rows))
                                     do.call(rbind, full_stats_rows) else NULL
    attr(out, "n_per_condition") <- n_per_cond_df

    # ---- Verbose: estimated min fig.width --------------------------------
    if (isTRUE(verbose)) {
        axis_width_cm   <- 1.4   # accounts for the row-label y-axis title in col 1
        legend_width_cm <- if (legend_position == "right") 3.0 else 0
        gap_cm          <- panel_gap_pt * 0.03528
        col_width_cm    <- panel_width_cm + 2 * gap_cm
        total_cm        <- n_cols_grid * col_width_cm + axis_width_cm + legend_width_cm
        total_in        <- total_cm / 2.54
        message(sprintf(
            "plotPbExprsDiffStrip: estimated min fig.width = %.1f in (%.1f cm) [n_cols=%d, panel_width_cm=%g, legend='%s']. Increase if PDFs look clipped.",
            total_in, total_cm, n_cols_grid, panel_width_cm, legend_position
        ))
    }

    out
}

# Attach CATALYST's namespace so .check_sce / .check_k / .check_assay /
# .check_cd_factor / .get_features / .agg / cluster_ids resolve. Same
# pattern as plotPbExprsDiff and plotAbundancesDiffStrip.
