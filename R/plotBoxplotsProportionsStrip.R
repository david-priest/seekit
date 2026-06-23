# plotBoxplotsProportionsStrip
#
# Sibling of plotBoxplotsProportions (and plotBoxplotsProportions2) with the
# TRUE separated stats-strip layout: each (master_cluster, daughter_cluster)
# cell of the 2D facet grid is composed of two stacked sub-plots,
#
#   [ stats strip ]   <- brackets only, no axes, fixed-height across the grid
#   [   data plot ]   <- daughter-within-master proportion boxplots
#
# Modeled directly on plotPbExprsDiffStrip.R. The structural parallel:
#   plotPbExprsDiff:    cluster_id     x antigen             (with `value` = median expr)
#   plotBoxplotsProps:  master_cluster x daughter_cluster    (with `proportion` = % within master)
#
# Grid orientation:
#   swap = FALSE -> rows = master_cluster, cols = daughter_cluster
#   swap = TRUE  -> rows = daughter_cluster, cols = master_cluster
# (matches plotBoxplotsProportions / plotBoxplotsProportions2 convention)
#
# Grid-label handling (same idea as plotPbExprsDiffStrip):
#   - top-row cells own the column variable as a strip title;
#   - leftmost-column cells embed the row variable into the y-axis title;
#   - interior cells suppress both labels.
#
# Stats provenance:
#   - The full `external_stats` (un-filtered by hide_ns), restricted to plotted
#     cells, is stashed on the returned patchwork as attr(p, "source_stats")
#     so f2() writes exact p-values for every comparison (including ns) to a
#     single "stats" sheet. attr(p, "source_data") is the long proportion table.
#
# Dependencies (loaded by SETUP 1 of the calling .Rmd):
#   SingleCellExperiment, ggplot2, ggbeeswarm, patchwork, ggh4x, scales, grid,
#   dplyr, tidyr, rlang. Optional: ggtext (when use_richtext = TRUE).

plotBoxplotsProportionsStrip <- function(
    sce,
    master_cluster_column,
    daughter_cluster_column,
    sample_id_column,
    group_column,
    meta              = c(sample_id_column, group_column),
    master_order      = NULL,
    daughter_order    = NULL,
    group_order       = NULL,
    excluded_masters  = NULL,
    excluded_daughters = NULL,
    fill_palette      = NULL,
    point_size        = 2,
    facet_ratio       = 1,
    panel_width_cm    = 3,
    swap              = FALSE,
    # Stats
    external_stats    = NULL,
    show_stats        = TRUE,
    hide_ns           = TRUE,
    # Multiple-testing correction view. Three options:
    #   "per_pair"    -> diffcyt default; BH across clusters within each pair
    #                   (m = clusters in that pair-call). Uses p.adj.signif.
    #   "per_cluster" -> recommended (per PI); BH across pairs within each
    #                   (master, daughter) only. Doesn't correct across other
    #                   masters or daughters. Uses p.adj_per_cluster.signif.
    #   "global"      -> BH across the full master × daughter × pairs family.
    #                   Uses p.adj_global.signif.
    # Convenience switch over `label_col`; if `label_col` is set explicitly
    # to anything other than the default, it wins.
    correction        = c("per_pair", "per_cluster", "global"),
    # Asterisk / label styling
    label_col         = "p.adj.signif",
    label_size        = 3,
    label_fontface    = "plain",
    use_richtext      = TRUE,
    asterisk_pt_multiplier = 1.4,
    asterisk_vjust    = 0.3,
    label_vjust       = 0,
    asterisk_y_offset = 0,
    bracket_color     = "black",
    bracket_size      = 0.4,
    bracket_tip_npc   = 0.20,
    label_nudge_npc   = 0.05,
    # Strip layout
    strip_height_fraction = 0.25,
    panel_gap_pt          = 8,
    # Y-axis
    y_expand_low_mult  = 0.02,
    y_expand_high_mult = 0.05,
    tight_top_axis           = TRUE,
    tight_top_axis_overhang  = 0.02,
    # Title placement / typography
    title_above_strip = TRUE,
    title_size        = 10,
    axis_text_size    = 8,
    row_label_size    = NULL,
    # Y-axis label
    y_axis_root       = "Proportion [%]",
    # Theme
    nature_style      = TRUE,
    legend_position   = "right",
    # Diagnostic
    verbose           = TRUE
) {
    # ---- Package dependencies --------------------------------------------
    if (!requireNamespace("patchwork", quietly = TRUE))
        stop("patchwork is required: install.packages('patchwork').")
    if (!requireNamespace("ggh4x", quietly = TRUE))
        stop("ggh4x is required (force_panelsizes): install.packages('ggh4x').")
    if (!requireNamespace("tidyr", quietly = TRUE))
        stop("tidyr is required (complete): install.packages('tidyr').")
    if (isTRUE(use_richtext) && !requireNamespace("ggtext", quietly = TRUE))
        stop("ggtext is required when use_richtext = TRUE: install.packages('ggtext').")
    if (is.null(row_label_size)) row_label_size <- title_size

    # Correction-view convenience switch (see param doc above).
    correction <- match.arg(correction)
    if (identical(label_col, "p.adj.signif")) {
        label_col <- switch(correction,
                            per_pair    = "p.adj.signif",
                            per_cluster = "p.adj_per_cluster.signif",
                            global      = "p.adj_global.signif")
    }

    # ---- Pull colData ----------------------------------------------------
    df <- as.data.frame(SingleCellExperiment::colData(sce))

    # ---- Factor levels (matches plotBoxplotsProportions[2] conventions) ----
    df[[master_cluster_column]]   <- factor(df[[master_cluster_column]],
                                            levels = if (!is.null(master_order)) master_order
                                                     else unique(df[[master_cluster_column]]))
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]],
                                            levels = if (!is.null(daughter_order)) daughter_order
                                                     else unique(df[[daughter_cluster_column]]))
    df[[sample_id_column]]        <- factor(df[[sample_id_column]],
                                            levels = unique(df[[sample_id_column]]))
    df[[group_column]]            <- factor(df[[group_column]],
                                            levels = if (!is.null(group_order)) group_order
                                                     else unique(df[[group_column]]))

    if (!is.null(excluded_masters)) {
        df <- df[!as.character(df[[master_cluster_column]]) %in% excluded_masters, , drop = FALSE]
        df[[master_cluster_column]] <- droplevels(df[[master_cluster_column]])
    }
    if (!is.null(excluded_daughters)) {
        df <- df[!as.character(df[[daughter_cluster_column]]) %in% excluded_daughters, , drop = FALSE]
        df[[daughter_cluster_column]] <- droplevels(df[[daughter_cluster_column]])
    }

    masters   <- levels(df[[master_cluster_column]])
    daughters <- levels(df[[daughter_cluster_column]])

    # ---- Proportion computation (daughter within master per sample) ------
    count_df <- df %>%
        dplyr::group_by_at(dplyr::vars(dplyr::all_of(
            c(sample_id_column, master_cluster_column, daughter_cluster_column, group_column)
        ))) %>%
        dplyr::summarise(count = dplyr::n(), .groups = "drop")

    complete_df <- count_df %>%
        tidyr::complete(
            tidyr::nesting(!!rlang::sym(sample_id_column),
                           !!rlang::sym(master_cluster_column),
                           !!rlang::sym(group_column)),
            !!rlang::sym(daughter_cluster_column),
            fill = list(count = 0)
        )

    proportion_df <- complete_df %>%
        dplyr::group_by_at(dplyr::vars(dplyr::all_of(
            c(sample_id_column, master_cluster_column, group_column)
        ))) %>%
        dplyr::mutate(proportion = count / sum(count) * 100) %>%
        dplyr::ungroup() %>%
        as.data.frame()

    proportion_df[[master_cluster_column]]   <- factor(proportion_df[[master_cluster_column]],   levels = masters)
    proportion_df[[daughter_cluster_column]] <- factor(proportion_df[[daughter_cluster_column]], levels = daughters)
    proportion_df[[group_column]]            <- factor(proportion_df[[group_column]],            levels = levels(df[[group_column]]))
    group_levels <- levels(proportion_df[[group_column]])

    # ---- n samples per condition (Immunity figure legend reporting) -------
    n_per_cond_df <- proportion_df %>%
        dplyr::distinct(!!rlang::sym(sample_id_column), !!rlang::sym(group_column)) %>%
        dplyr::count(!!rlang::sym(group_column), name = "n_samples") %>%
        as.data.frame()
    message("plotBoxplotsProportionsStrip: n samples per ", group_column, ":")
    print(n_per_cond_df, row.names = FALSE)

    # ---- Per-cell stats + global max_brackets ----------------------------
    cell_stats <- list()
    full_stats <- NULL
    if (isTRUE(show_stats) && !is.null(external_stats)) {
        es <- as.data.frame(external_stats)
        keep <- as.character(es[[master_cluster_column]])   %in% masters &
                as.character(es[[daughter_cluster_column]]) %in% daughters
        full_stats <- es[keep, , drop = FALSE]
        for (mc in masters) for (dc in daughters) {
            sub <- full_stats[as.character(full_stats[[master_cluster_column]])   == mc &
                              as.character(full_stats[[daughter_cluster_column]]) == dc, , drop = FALSE]
            if (isTRUE(hide_ns)) {
                sub <- sub[as.character(sub[[label_col]]) != "ns", , drop = FALSE]
            }
            cell_stats[[paste(mc, dc, sep = "\037")]] <- sub
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

    # ---- Label transform (asterisk styling, same as siblings) -----------
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

    # ---- Per-cell builders -----------------------------------------------
    build_data_plot <- function(df_cell, row_lbl, col_lbl,
                                show_row_label, show_col_title) {
        y_lab <- if (isTRUE(show_row_label)) paste0(row_lbl, "\n", y_axis_root) else NULL

        p <- ggplot2::ggplot(
                df_cell,
                ggplot2::aes(x    = .data[[group_column]],
                             y    = proportion,
                             fill = .data[[group_column]])
            ) +
            ggplot2::geom_boxplot(
                color = "grey16", linewidth = 0.5, width = 0.75,
                alpha = 0.8, outlier.color = NA,
                show.legend = TRUE
            ) +
            ggbeeswarm::geom_quasirandom(
                shape = 21, fill = "grey84", size = point_size,
                width = 0.2, alpha = 0.8
            ) +
            ggplot2::labs(
                x = NULL,
                y = y_lab,
                title = if (isTRUE(show_col_title) && !isTRUE(title_above_strip))
                            col_lbl else NULL
            ) +
            (if (isTRUE(tight_top_axis)) {
                y_max <- max(df_cell$proportion, na.rm = TRUE)
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
            p <- p + ggplot2::scale_fill_manual(values = fill_palette, name = group_column, drop = FALSE)
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
                    label.size = NA, fill = NA,
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

    # ---- Build 2D grid of composites -------------------------------------
    # swap = FALSE -> rows = master, cols = daughter (matches the upstream
    #                 plotBoxplotsProportions default).
    # swap = TRUE  -> rows = daughter, cols = master  (matches the upstream
    #                 swap = TRUE usage, e.g. isotype proportions chunk).
    if (!isTRUE(swap)) {
        row_vals <- masters;   col_vals <- daughters
    } else {
        row_vals <- daughters; col_vals <- masters
    }
    n_cols_grid <- length(col_vals)

    composites <- list()
    idx <- 0L
    for (rv in row_vals) for (cv in col_vals) {
        idx <- idx + 1L
        if (!isTRUE(swap)) { mc <- rv; dc <- cv } else { mc <- cv; dc <- rv }

        df_cell <- proportion_df[as.character(proportion_df[[master_cluster_column]])   == mc &
                                 as.character(proportion_df[[daughter_cluster_column]]) == dc, , drop = FALSE]
        stat_cell <- cell_stats[[paste(mc, dc, sep = "\037")]]

        row_pos <- ((idx - 1L) %/% n_cols_grid) + 1L
        col_pos <- ((idx - 1L) %%  n_cols_grid) + 1L
        show_col_title <- (row_pos == 1L)
        show_row_label <- (col_pos == 1L)

        if (nrow(df_cell) == 0) {
            df_cell <- data.frame(proportion = c(0, 1))
            df_cell[[group_column]] <- factor(group_levels[c(1, 1)], levels = group_levels)
        }

        p_data  <- build_data_plot(df_cell, row_lbl = rv, col_lbl = cv,
                                   show_row_label = show_row_label,
                                   show_col_title = show_col_title)
        p_strip <- build_strip_plot(stat_cell, col_lbl = cv,
                                    show_col_title = show_col_title)

        composites[[idx]] <- p_strip / p_data +
            patchwork::plot_layout(heights = c(strip_h, data_h))
    }

    gap <- panel_gap_pt
    out <- patchwork::wrap_plots(composites, ncol = n_cols_grid) &
        ggplot2::theme(
            legend.position = legend_position,
            plot.margin     = ggplot2::margin(t = gap, r = gap, b = gap, l = gap, unit = "pt")
        )
    if (legend_position != "none") {
        out <- out + patchwork::plot_layout(guides = "collect")
    }

    # ---- Stash long-format data + full unfiltered stats ------------------
    attr(out, "source_data")  <- proportion_df
    attr(out, "n_per_condition") <- n_per_cond_df
    attr(out, "source_stats") <- full_stats

    if (isTRUE(verbose)) {
        axis_width_cm   <- 1.4
        legend_width_cm <- if (legend_position == "right") 3.0 else 0
        gap_cm          <- panel_gap_pt * 0.03528
        col_width_cm    <- panel_width_cm + 2 * gap_cm
        total_cm        <- n_cols_grid * col_width_cm + axis_width_cm + legend_width_cm
        total_in        <- total_cm / 2.54
        message(sprintf(
            "plotBoxplotsProportionsStrip: estimated min fig.width = %.1f in (%.1f cm) [n_cols=%d, panel_width_cm=%g, legend='%s']. Increase if PDFs look clipped.",
            total_in, total_cm, n_cols_grid, panel_width_cm, legend_position
        ))
    }

    out
}
