# plotBoxplotsRatio
#
# Per-master-cluster boxplots of the ratio of two daughter-cluster groups.
# Cousin of plotBoxplotsProportions (which plots daughter-proportion boxplots
# faceted by master x daughter) and plotCD4CD8Ratio (which plots a single
# ratio without per-master-cluster faceting).
#
# Use case: for each B-cell merging8_B cluster (Activated / Atypical / ...),
# plot the per-sample ratio of IgM+ to IgG+ cells (the "daughter clusters",
# living in the `isotype` colData column).
#
# - numerator_cluster / denominator_cluster can each be a single cluster label
#   or a character vector (multiple labels get summed before division).
# - Adds a dotted unity reference line at y = 1.
# - Supports log-scale y axis (with miny pseudocount shift, same convention
#   as plotAbundancesDiff / plotCD4CD8Ratio).

plotBoxplotsRatio <- function(
    sce,
    master_cluster_column,
    daughter_cluster_column,
    numerator_cluster,
    denominator_cluster,
    group_column,
    sample_id_column,
    meta            = NULL,            # additional colData columns to preserve per sample (e.g. c("sample_id","patient_id","condition3"))
    group_order     = NULL,
    group_palette   = NULL,
    facet_ratio     = 1.5,
    panel_spacing   = 1,
    point_size      = 1.5,
    n_cols          = NULL,
    pseudocount     = 0,
    show_stats      = FALSE,
    step_increase   = 0,
    external_stats  = NULL,
    log             = FALSE,
    miny            = 0.01,
    maxy            = NA,
    hide_ns         = FALSE,      # TRUE hides non-significant brackets
    # ---- NPC (Prism-style brackets above panels) ----
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

    # ---- Input validation ------------------------------------------------
    required_cols <- c(master_cluster_column, daughter_cluster_column,
                       group_column, sample_id_column)
    missing_cols <- setdiff(required_cols, colnames(df))
    if (length(missing_cols) > 0) {
        stop("Missing required colData columns: ",
             paste(missing_cols, collapse = ", "))
    }
    if (!is.null(meta)) {
        missing_meta <- setdiff(meta, colnames(df))
        if (length(missing_meta) > 0) {
            stop("Missing meta columns in colData: ",
                 paste(missing_meta, collapse = ", "))
        }
    }
    if (!is.character(numerator_cluster) || length(numerator_cluster) < 1) {
        stop("numerator_cluster must be a non-empty character vector.")
    }
    if (!is.character(denominator_cluster) || length(denominator_cluster) < 1) {
        stop("denominator_cluster must be a non-empty character vector.")
    }
    if (length(intersect(numerator_cluster, denominator_cluster)) > 0) {
        stop("numerator_cluster and denominator_cluster must not overlap.")
    }
    if (!is.numeric(pseudocount) || length(pseudocount) != 1 || pseudocount < 0) {
        stop("pseudocount must be a single non-negative numeric value.")
    }

    # Coerce to character for stable matching
    df[[master_cluster_column]]   <- as.character(df[[master_cluster_column]])
    df[[daughter_cluster_column]] <- as.character(df[[daughter_cluster_column]])

    observed_daughters <- unique(df[[daughter_cluster_column]])
    missing_daughters  <- setdiff(c(numerator_cluster, denominator_cluster),
                                  observed_daughters)
    if (length(missing_daughters) > 0) {
        warning("daughter cluster(s) not found in '", daughter_cluster_column,
                "': ", paste(missing_daughters, collapse = ", "))
    }

    # ---- Tag each cell as numerator / denominator / other -----------------
    df$.daughter_role <- ifelse(
        df[[daughter_cluster_column]] %in% numerator_cluster,   "numerator",
        ifelse(df[[daughter_cluster_column]] %in% denominator_cluster,
               "denominator", NA_character_)
    )
    df <- df[!is.na(df$.daughter_role), , drop = FALSE]
    if (nrow(df) == 0) {
        stop("No cells matched numerator_cluster or denominator_cluster.")
    }

    # ---- Per-(sample, master, role) counts and ratio ----------------------
    # `meta` columns are sample-level (constant within sample_id) and get added
    # to the groupby so they propagate through to the ratio_df. Duplicates with
    # the core group columns are collapsed.
    group_cols <- unique(c(sample_id_column, master_cluster_column,
                           group_column, meta))
    count_df <- df %>%
        dplyr::group_by(dplyr::across(
            dplyr::all_of(c(group_cols, ".daughter_role"))
        )) %>%
        dplyr::summarise(count = dplyr::n(), .groups = "drop")

    # Make sure every (sample, master) has both roles (fill zeros)
    count_df <- count_df %>%
        tidyr::complete(
            tidyr::nesting(!!!rlang::syms(group_cols)),
            .daughter_role = c("numerator", "denominator"),
            fill = list(count = 0)
        )

    ratio_df <- count_df %>%
        tidyr::pivot_wider(names_from  = .daughter_role,
                           values_from = count,
                           values_fill = 0)

    if (!"numerator"   %in% colnames(ratio_df)) ratio_df$numerator   <- 0
    if (!"denominator" %in% colnames(ratio_df)) ratio_df$denominator <- 0

    if (pseudocount == 0) {
        ratio_df$ratio <- ifelse(ratio_df$denominator == 0, NA_real_,
                                 ratio_df$numerator / ratio_df$denominator)
    } else {
        ratio_df$ratio <- (ratio_df$numerator   + pseudocount) /
                          (ratio_df$denominator + pseudocount)
    }

    dropped_n <- sum(is.na(ratio_df$ratio))
    if (dropped_n > 0) {
        warning(dropped_n, " sample/master combinations have zero ",
                "denominator cells; ratio set to NA for those rows.")
    }
    ratio_df <- ratio_df %>% dplyr::filter(!is.na(ratio))

    # ---- Factor levels ----------------------------------------------------
    if (!is.null(group_order)) {
        ratio_df[[group_column]] <- factor(ratio_df[[group_column]],
                                           levels = group_order)
    }

    if (is.null(group_palette)) {
        group_palette <- RColorBrewer::brewer.pal(
            min(length(unique(ratio_df[[group_column]])), 12), "Paired"
        )
    }

    # ---- Log-scale shift (matches plotAbundancesDiff convention) ---------
    if (log) {
        ratio_df$ratio <- ratio_df$ratio + miny
    }

    ratio_out <<- ratio_df

    # ---- Axis labels ------------------------------------------------------
    numerator_label   <- paste(numerator_cluster,   collapse = "+")
    denominator_label <- paste(denominator_cluster, collapse = "+")
    y_label <- paste0(numerator_label, " / ", denominator_label, " ratio")

    y_breaks_fun <- function(x) {
        sort(unique(c(scales::breaks_pretty(n = 5)(x), 1)))
    }

    # ---- Plot -------------------------------------------------------------
    p <- ggplot2::ggplot(
            ratio_df,
            ggplot2::aes(x    = !!rlang::sym(group_column),
                         y    = ratio,
                         fill = !!rlang::sym(group_column))
        ) +
        ggplot2::labs(x = NULL, y = y_label) +
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
        # Unity reference line (ratio == 1). Drawn before boxes so it sits beneath.
        ggplot2::geom_hline(yintercept = 1, linetype = "dotted",
                            color = "grey40", linewidth = 0.5) +
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
        ggplot2::facet_wrap(
            ggplot2::vars(!!rlang::sym(master_cluster_column)),
            scales = "free_y",
            ncol   = n_cols
        )

    if (log) {
        p <- p +
            ggplot2::scale_y_continuous(
                trans  = "log10",
                limits = c(miny, maxy),
                breaks = c(0.01, 0.1, 1, 10, 100),
                labels = c(0.01, 0.1, 1, 10, 100)
            ) +
            ggplot2::annotation_logticks(base = 10, sides = "l", outside = TRUE) +
            ggplot2::coord_cartesian(clip = "off") +
            ggplot2::theme(axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)))
    } else {
        p <- p +
            ggplot2::scale_y_continuous(limits = c(0, NA), breaks = y_breaks_fun) +
            ggplot2::expand_limits(y = 1) +
            ggplot2::coord_cartesian(clip = "off")
    }

    # ---- Stats ------------------------------------------------------------
    if (show_stats) {
        # --- NPC (Prism-style) path: route brackets through add_pvalue_npc ---
        if (identical(stats_style, "npc")) {
            stat_input <- ratio_df %>% dplyr::filter(!is.na(ratio))

            # Build the stat.test for add_pvalue_npc. external_stats is rare for
            # this function (no diffcyt parallel), so we run dunn_test inline.
            if (is.null(external_stats)) {
                # Per-master_cluster Dunn's test, with safety try/catch (some
                # master clusters may have too few samples or only one group).
                split_input <- split(stat_input, stat_input[[master_cluster_column]])
                stat_chunks <- lapply(names(split_input), function(grp) {
                    d <- split_input[[grp]]
                    if (dplyr::n_distinct(d[[group_column]]) < 2 || nrow(d) < 3) {
                        return(NULL)
                    }
                    res <- tryCatch(
                        rstatix::dunn_test(d,
                                           as.formula(paste("ratio ~", group_column)),
                                           p.adjust.method = "holm"),
                        error = function(e) NULL
                    )
                    if (is.null(res) || nrow(res) == 0) return(NULL)
                    res[[master_cluster_column]] <- grp
                    res
                })
                stat_chunks <- stat_chunks[!vapply(stat_chunks, is.null, logical(1))]
                if (length(stat_chunks) == 0) {
                    warning("No master_cluster had enough data for stats; skipping bracket layer.")
                    npc_stat <- NULL
                } else {
                    npc_stat <- dplyr::bind_rows(stat_chunks)
                }
            } else {
                npc_stat <- external_stats
            }

            # Per-master_cluster max ratio for panel_max
            panel_max_df <- stat_input %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(master_cluster_column))) %>%
                dplyr::summarise(max_val = max(ratio, na.rm = TRUE), .groups = "drop") %>%
                dplyr::mutate(max_val = ifelse(max_val == 0 | is.na(max_val), 1, max_val))

            grp_levels <- if (!is.null(group_order)) {
                group_order
            } else if (is.factor(ratio_df[[group_column]])) {
                levels(droplevels(ratio_df[[group_column]]))
            } else {
                sort(unique(as.character(ratio_df[[group_column]])))
            }

            statout <<- npc_stat
            if (!is.null(npc_stat) && nrow(npc_stat) > 0) {
                p <- p + add_pvalue_npc(
                    stat.test              = npc_stat,
                    panel_max              = panel_max_df,
                    facet_var              = master_cluster_column,
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
        stat_input <- ratio_df %>% dplyr::filter(!is.na(ratio))

        # Drop master_cluster groups that don't have enough data for stats:
        #   - need >= 2 distinct levels of group_column to make ANY pair
        #   - need >= 3 samples total (rstatix dunn_test needs at least that)
        # Otherwise rstatix throws "New column has 6 rows. .data has 0 rows".
        master_summary <- stat_input %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(master_cluster_column))) %>%
            dplyr::summarise(
                n_groups  = dplyr::n_distinct(!!rlang::sym(group_column)),
                n_samples = dplyr::n(),
                .groups   = "drop"
            )
        valid_masters <- master_summary %>%
            dplyr::filter(n_groups >= 2, n_samples >= 3) %>%
            dplyr::pull(!!rlang::sym(master_cluster_column))
        dropped_masters <- setdiff(master_summary[[master_cluster_column]], valid_masters)
        if (length(dropped_masters) > 0) {
            warning("Skipping stats for master_cluster group(s) with insufficient data: ",
                    paste(dropped_masters, collapse = ", "))
        }
        stat_input <- stat_input %>%
            dplyr::filter(!!rlang::sym(master_cluster_column) %in% valid_masters)

        if (nrow(stat_input) == 0) {
            warning("No master_cluster groups had enough data for stats; skipping stat layer.")
            stat.test <- NULL
        } else if (!is.null(external_stats)) {
            external_stats[[master_cluster_column]] <- as.character(
                external_stats[[master_cluster_column]]
            )
            # Keep only external stats for valid masters too
            external_stats <- external_stats %>%
                dplyr::filter(!!rlang::sym(master_cluster_column) %in% valid_masters)

            eff_step <- if (step_increase == 0) 0.1 else step_increase

            max_y_df <- stat_input %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(master_cluster_column))) %>%
                dplyr::summarise(max_val = max(ratio, na.rm = TRUE), .groups = "drop") %>%
                dplyr::mutate(max_val = ifelse(max_val == 0 | is.na(max_val), 1, max_val))

            stat.test <- external_stats %>%
                dplyr::left_join(max_y_df, by = master_cluster_column) %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(master_cluster_column))) %>%
                dplyr::mutate(
                    y.position = max_val * 1.1 + (dplyr::row_number() - 1) * (max_val * eff_step)
                ) %>%
                dplyr::ungroup()
        } else {
            # Wrap each per-group dunn_test in tryCatch so a single failing group
            # doesn't blow up the whole stat computation.
            split_input <- split(stat_input, stat_input[[master_cluster_column]])
            per_group_stats <- lapply(names(split_input), function(grp) {
                d <- split_input[[grp]]
                res <- tryCatch(
                    rstatix::dunn_test(d,
                                       as.formula(paste("ratio ~", group_column)),
                                       p.adjust.method = "holm"),
                    error = function(e) NULL
                )
                if (is.null(res) || nrow(res) == 0) return(NULL)
                res <- tryCatch(
                    rstatix::add_y_position(res,
                                            step.increase = step_increase,
                                            fun = "max",
                                            data = d,
                                            formula = as.formula(paste("ratio ~", group_column))),
                    error = function(e) NULL
                )
                if (is.null(res) || nrow(res) == 0) return(NULL)
                res[[master_cluster_column]] <- grp
                res
            })
            per_group_stats <- per_group_stats[!vapply(per_group_stats, is.null, logical(1))]
            if (length(per_group_stats) == 0) {
                warning("Stat computation failed for every master_cluster group; skipping stat layer.")
                stat.test <- NULL
            } else {
                stat.test <- dplyr::bind_rows(per_group_stats)
            }
        }

        statout <<- stat.test

        if (!is.null(stat.test) && nrow(stat.test) > 0) {
            # inherit.aes = FALSE + explicit group1/group2 mapping is needed because
            # the parent ggplot is keyed on x = <group_column> (e.g. condition3) but
            # stat.test only has group1/group2 columns. Without the override,
            # ggplot complains "object 'condition3' not found" when the stat
            # layer tries to inherit the plot's aes.
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
    }

    return(p)
}
