# plotBoxplotsProportions2
#
# Sibling of plotBoxplotsProportions with an added `log` option for plotting
# proportions on a log10 y-axis. Default behavior matches v1 (linear, no stats).
#
# Key differences from v1:
#   - `log = FALSE` (default) -> identical behavior to plotBoxplotsProportions.
#   - `log = TRUE` -> y-axis on log10 scale. Proportions are shifted by `miny`
#     before transform (additive pseudocount, same convention as
#     plotAbundancesDiff / plotCD4CD8Ratio), so zero counts plot at miny instead
#     of -Inf.
#   - `show_stats = FALSE` is the new default. Stats brackets on log axes are
#     fiddly (the v1 y.position math is additive and lands the brackets in
#     wrong places under log10); you can still pass show_stats = TRUE but you
#     have been warned.
#
# Notes on what *didn't* work in earlier attempts and why this version avoids it:
#   - scale_y_continuous(limits = c(miny, maxy)) is a global hard limit and
#     fights with ggh4x::facet_grid2(independent = "y"). DROPPED -- let each
#     facet auto-scale within the log10 transform.
#   - annotation_logticks(outside = TRUE) requires coord_cartesian(clip = "off")
#     which conflicts with ggh4x's coord. DROPPED.
#   - coord_cartesian() overrides ggh4x's faceting coord. DROPPED.

plotBoxplotsProportions2 <- function(
    sce,
    master_cluster_column,
    daughter_cluster_column,
    meta,
    group_column,
    sample_id_column,
    plot_order_df  = NULL,
    daughter_order = NULL,
    group_order    = NULL,
    group_palette  = NULL,
    facet_ratio    = 1.5,
    panel_spacing  = 1,
    point_size     = 1.5,
    show_stats     = FALSE,        # off by default; stats can wreck log axes
    step_increase  = 0,
    external_stats = NULL,
    swap           = FALSE,
    log            = FALSE,        # NEW
    miny           = 0.01,         # NEW: additive pseudocount for log path
    maxy           = NA,           # NEW: reserved, currently unused
    hide_ns        = FALSE,        # TRUE hides non-significant brackets
    # ---- NPC (Prism-style brackets, 2D facet) ----
    stats_style              = c("data", "npc"),
    npc_bracket_top_offset   = 0.25,
    npc_bracket_step         = 0.20,
    npc_label_size           = 2.5,
    npc_richtext             = TRUE,
    npc_asterisk_pt_multiplier = 1.4,
    npc_panel_spacing_y_lines  = 3,
    npc_plot_margin_top_pt     = 30,
    npc_nature_style           = FALSE
) {
    stats_style <- match.arg(stats_style)
    df <- as.data.frame(SingleCellExperiment::colData(sce))

    # Factor levels (matches v1)
    df[[master_cluster_column]]   <- factor(df[[master_cluster_column]],
                                            levels = unique(df[[master_cluster_column]]))
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]],
                                            levels = unique(df[[daughter_cluster_column]]))
    df[[sample_id_column]]        <- factor(df[[sample_id_column]],
                                            levels = unique(df[[sample_id_column]]))
    df[[group_column]]            <- factor(df[[group_column]],
                                            levels = if (!is.null(group_order)) group_order
                                                     else unique(df[[group_column]]))
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]],
                                            levels = if (!is.null(daughter_order)) daughter_order
                                                     else unique(df[[daughter_cluster_column]]))

    if (is.null(group_palette))
        group_palette <- RColorBrewer::brewer.pal(min(length(levels(df[[group_column]])), 12), "Paired")

    # Per-sample-master-daughter counts -> per-sample-master proportions
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
        dplyr::ungroup()

    summary_df <- proportion_df %>% dplyr::rename(cell_count = count)

    # ---- LOG: additive pseudocount shift before any transform ----
    # Same convention as plotAbundancesDiff(): zeros land at miny instead of -Inf.
    if (log) {
        summary_df$proportion <- summary_df$proportion + miny
    }

    summary_out <<- summary_df

    if (log && show_stats) {
        warning("plotBoxplotsProportions2(): show_stats = TRUE with log = TRUE -- ",
                "bracket positions are computed additively and may land in ",
                "unexpected places on the log axis. Default show_stats = FALSE ",
                "for log plots is recommended.")
    }

    # ---- Base plot ----
    p <- ggplot2::ggplot(summary_df,
                         ggplot2::aes(x    = !!rlang::sym(group_column),
                                      y    = proportion,
                                      fill = !!rlang::sym(group_column))) +
        ggplot2::geom_boxplot(
            color        = "grey16",
            position     = ggplot2::position_dodge(),
            size         = 0.5,
            alpha        = 0.8,
            outlier.color = NA,
            show.legend  = TRUE
        ) +
        ggbeeswarm::geom_quasirandom(
            fill  = "grey84",
            size  = point_size,
            width = 0.2,
            shape = 21,
            alpha = 0.8
        ) +
        ggplot2::scale_fill_manual(values = group_palette, name = group_column) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid        = ggplot2::element_blank(),
            strip.text        = ggplot2::element_text(face = "bold"),
            strip.background  = ggplot2::element_rect(fill = NA, color = NA),
            axis.text         = ggplot2::element_text(color = "black"),
            aspect.ratio      = facet_ratio,
            axis.text.x       = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
            axis.ticks        = ggplot2::element_line(color = "black"),
            panel.border      = ggplot2::element_rect(color = "black", fill = NA, size = 0.5),
            panel.spacing     = ggplot2::unit(panel_spacing, "lines"),
            legend.key.height = ggplot2::unit(1.5, "lines")
        ) +
        ggplot2::labs(
            x = NULL,
            y = if (log) "Proportion [%] (log10)" else "Proportion [%]"
        )

    # ---- Y-axis transform ----
    # Use ONLY the trans, no hard limits, no coord override. ggh4x's
    # independent="y" then auto-fits each facet to its own log range.
    if (log) {
        p <- p + ggplot2::scale_y_continuous(trans = "log10")
    }

    # ---- Faceting ----
    if (swap) {
        p <- p + ggh4x::facet_grid2(
            cols          = ggplot2::vars(!!rlang::sym(master_cluster_column)),
            rows          = ggplot2::vars(!!rlang::sym(daughter_cluster_column)),
            scales        = "free",
            independent   = "y",
            remove_labels = "x"
        )
    } else {
        p <- p + ggh4x::facet_grid2(
            rows          = ggplot2::vars(!!rlang::sym(master_cluster_column)),
            cols          = ggplot2::vars(!!rlang::sym(daughter_cluster_column)),
            scales        = "free",
            independent   = "y",
            remove_labels = "x"
        )
    }

    # ---- Stats (default off) ----
    if (show_stats) {
        # --- NPC (Prism-style) path: 2D facet via composite key ---
        if (identical(stats_style, "npc")) {
            panel_max_df <- summary_df %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(c(master_cluster_column, daughter_cluster_column)))) %>%
                dplyr::summarise(max_val = max(proportion, na.rm = TRUE), .groups = "drop") %>%
                dplyr::mutate(max_val = ifelse(max_val == 0 | is.na(max_val), 1, max_val))

            if (is.null(external_stats)) {
                warning("plotBoxplotsProportions2: stats_style='npc' currently requires external_stats. Skipping bracket layer.")
                npc_stat <- NULL
            } else {
                npc_stat <- external_stats
            }

            grp_levels <- if (!is.null(group_order)) {
                group_order
            } else if (is.factor(df[[group_column]])) {
                levels(droplevels(df[[group_column]]))
            } else {
                sort(unique(as.character(df[[group_column]])))
            }

            statout <<- npc_stat
            if (!is.null(npc_stat) && nrow(npc_stat) > 0) {
                p <- p + add_pvalue_npc(
                    stat.test              = npc_stat,
                    panel_max              = panel_max_df,
                    facet_var              = c(master_cluster_column, daughter_cluster_column),
                    max_col                = "max_val",
                    group_levels           = grp_levels,
                    hide_ns                = hide_ns,
                    bracket_top_offset     = npc_bracket_top_offset,
                    bracket_step           = npc_bracket_step,
                    label_size             = npc_label_size,
                    use_richtext           = npc_richtext,
                    asterisk_pt_multiplier = npc_asterisk_pt_multiplier,
                    panel_spacing_y_lines  = npc_panel_spacing_y_lines,
                    plot_margin_top_pt     = npc_plot_margin_top_pt,
                    nature_style           = npc_nature_style
                )
            }
            return(p)
        }

        # --- Legacy data-space stats path ---
        if (!is.null(external_stats)) {
            external_stats$isotype <- as.factor(as.character(external_stats$isotype))
            summary_df$isotype     <- as.factor(as.character(summary_df$isotype))

            external_stats[[master_cluster_column]] <- factor(
                external_stats[[master_cluster_column]],
                levels = levels(df[[master_cluster_column]])
            )
            summary_df[[master_cluster_column]] <- factor(
                summary_df[[master_cluster_column]],
                levels = levels(df[[master_cluster_column]])
            )

            if (any(is.na(summary_df$proportion))) {
                stop("Missing values found in 'proportion' column of summary_df")
            }

            eff_step <- if (step_increase == 0) 0.1 else step_increase

            max_y_df <- summary_df %>%
                dplyr::group_by_at(dplyr::vars(dplyr::all_of(
                    c(master_cluster_column, daughter_cluster_column)
                ))) %>%
                dplyr::summarise(max_val = max(proportion, na.rm = TRUE), .groups = "drop") %>%
                dplyr::mutate(max_val = ifelse(max_val == 0 | is.na(max_val), 1, max_val))

            stat.test <- external_stats %>%
                dplyr::left_join(max_y_df,
                                 by = c(master_cluster_column, daughter_cluster_column)) %>%
                dplyr::group_by_at(dplyr::vars(dplyr::all_of(
                    c(master_cluster_column, daughter_cluster_column)
                ))) %>%
                dplyr::mutate(
                    # Multiplicative spacing on log; additive on linear.
                    y.position = if (log) {
                        max_val * (1.5 ^ (1 + (dplyr::row_number() - 1) * eff_step))
                    } else {
                        max_val * 1.1 + (dplyr::row_number() - 1) * (max_val * eff_step)
                    }
                ) %>%
                dplyr::ungroup()
        } else {
            stat.test <- summary_df %>%
                dplyr::group_by_at(dplyr::vars(dplyr::all_of(
                    c(master_cluster_column, daughter_cluster_column)
                ))) %>%
                rstatix::dunn_test(
                    formula = as.formula(paste("proportion ~", group_column)),
                    p.adjust.method = "holm"
                ) %>%
                rstatix::add_y_position(scales = "free", step.increase = step_increase)
        }

        statout <<- stat.test

        p <- p + ggpubr::stat_pvalue_manual(
            stat.test,
            label       = "p.adj.signif",
            tip.length  = 0.01,
            hide.ns     = hide_ns,
            size        = 5,
            y.position  = "y.position",
            mapping     = ggplot2::aes(x = group1, xend = group2),
            inherit.aes = FALSE
        )
    }

    return(p)
}
