library(SingleCellExperiment)
library(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(tidyr)
library(rstatix)

plotCD4CD8Ratio <- function(
    x,
    k = "merging1",
    group_by = "condition3",
    patient_by = "sample_id",
    shape_by = NULL,
    numerator_cluster = "CD4",
    denominator_cluster = "CD8",
    cd4_label = NULL,
    cd8_label = NULL,
    pseudocount = 0,
    point_size = 2,
    panel_spacing = 1,
    facet_ratio = 1.5,
    log = FALSE,
    miny = 0.01,
    maxy = NA,
    step_increase = 0,
    external_stats = NULL,
    show_stats = FALSE,
    hide_ns = FALSE,         # TRUE hides non-significant brackets
    # ---- NPC (Prism-style brackets, single-panel) ----
    stats_style              = c("data", "npc"),
    npc_bracket_top_offset   = 0.25,
    npc_bracket_step         = 0.20,
    npc_label_size           = 4,
    npc_richtext             = TRUE,
    npc_asterisk_pt_multiplier = 1.4,
    npc_panel_spacing_y_lines  = 3,
    npc_plot_margin_top_pt     = 30,
    npc_nature_style           = FALSE
) {
    stats_style <- match.arg(stats_style)
    if (!inherits(x, "SingleCellExperiment")) {
        stop("x must be a SingleCellExperiment object.")
    }

    df <- as.data.frame(SingleCellExperiment::colData(x))

    if (!(k %in% colnames(df))) {
        resolved_k <- NULL

        # Preferred route: derive k-level labels via the CATALYST-free resolver
        # (.wl_cluster_ids reads metadata(x)$cluster_codes; no CATALYST needed).
        resolved_k <- tryCatch(
            as.character(.wl_cluster_ids(x, k = k)),
            error = function(e) NULL
        )

        # Fallback route: map from metadata$cluster_codes using cluster_id (or row names).
        if (is.null(resolved_k)) {
            cluster_codes <- SingleCellExperiment::metadata(x)$cluster_codes

            if (!is.null(cluster_codes)) {
                if (!is.data.frame(cluster_codes)) {
                    cluster_codes <- as.data.frame(cluster_codes)
                }

                if (k %in% colnames(cluster_codes)) {
                    if ("cluster_id" %in% colnames(df)) {
                        id_vec <- as.character(df$cluster_id)
                        candidate_cols <- setdiff(colnames(cluster_codes), k)

                        key_col <- NULL
                        if (length(candidate_cols) > 0) {
                            key_hits <- vapply(candidate_cols, function(cc) {
                                vals <- as.character(cluster_codes[[cc]])
                                all(unique(id_vec) %in% vals)
                            }, logical(1))

                            if (any(key_hits)) {
                                key_col <- candidate_cols[which(key_hits)[1]]
                            }
                        }

                        if (!is.null(key_col)) {
                            idx <- match(id_vec, as.character(cluster_codes[[key_col]]))
                            resolved_k <- as.character(cluster_codes[[k]][idx])
                        } else if (!is.null(rownames(cluster_codes))) {
                            idx <- match(id_vec, rownames(cluster_codes))
                            resolved_k <- as.character(cluster_codes[[k]][idx])
                        }
                    }
                }
            }
        }

        if (is.null(resolved_k) || all(is.na(resolved_k))) {
            stop(
                "Missing required colData column: ", k,
                ". Could not derive it from .wl_cluster_ids or metadata(x)$cluster_codes."
            )
        }

        df[[k]] <- resolved_k
    }

    needed_cols <- c(k, group_by, patient_by)
    missing_cols <- setdiff(needed_cols, colnames(df))
    if (length(missing_cols) > 0) {
        stop("Missing required colData columns: ", paste(missing_cols, collapse = ", "))
    }

    if (!is.null(shape_by) && !(shape_by %in% colnames(df))) {
        stop("shape_by column not found in colData: ", shape_by)
    }

    if (!is.numeric(pseudocount) || length(pseudocount) != 1 || pseudocount < 0) {
        stop("pseudocount must be a single non-negative numeric value.")
    }

    # Backwards compatibility with the old argument names.
    if (!is.null(cd4_label)) {
        numerator_cluster <- cd4_label
    }
    if (!is.null(cd8_label)) {
        denominator_cluster <- cd8_label
    }

    if (!is.character(numerator_cluster) || length(numerator_cluster) != 1 || !nzchar(numerator_cluster)) {
        stop("numerator_cluster must be a single non-empty character value.")
    }
    if (!is.character(denominator_cluster) || length(denominator_cluster) != 1 || !nzchar(denominator_cluster)) {
        stop("denominator_cluster must be a single non-empty character value.")
    }
    if (identical(numerator_cluster, denominator_cluster)) {
        stop("numerator_cluster and denominator_cluster must be different.")
    }

    df$cluster_lab <- as.character(df[[k]])
    observed_clusters <- unique(df$cluster_lab)
    if (!(numerator_cluster %in% observed_clusters)) {
        warning("numerator_cluster not found in ", k, ": ", numerator_cluster)
    }
    if (!(denominator_cluster %in% observed_clusters)) {
        warning("denominator_cluster not found in ", k, ": ", denominator_cluster)
    }

    df <- df %>% dplyr::filter(cluster_lab %in% c(numerator_cluster, denominator_cluster))

    if (nrow(df) == 0) {
        stop(
            "No cells found with labels ", numerator_cluster,
            " and/or ", denominator_cluster,
            " in ", k, "."
        )
    }

    group_cols <- c(patient_by, group_by)
    if (!is.null(shape_by)) {
        group_cols <- c(group_cols, shape_by)
    }

    count_df <- df %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "cluster_lab")))) %>%
        dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")

    ratio_df <- count_df %>%
        tidyr::pivot_wider(names_from = cluster_lab, values_from = n_cells, values_fill = 0)

    if (!(numerator_cluster %in% colnames(ratio_df))) {
        ratio_df[[numerator_cluster]] <- 0
    }
    if (!(denominator_cluster %in% colnames(ratio_df))) {
        ratio_df[[denominator_cluster]] <- 0
    }

    numerator_n <- ratio_df[[numerator_cluster]]
    denominator_n <- ratio_df[[denominator_cluster]]

    if (pseudocount == 0) {
        ratio_df$ratio <- ifelse(denominator_n == 0, NA_real_, numerator_n / denominator_n)
    } else {
        ratio_df$ratio <- (numerator_n + pseudocount) / (denominator_n + pseudocount)
    }

    names(ratio_df)[names(ratio_df) == numerator_cluster] <- "numerator_n"
    names(ratio_df)[names(ratio_df) == denominator_cluster] <- "denominator_n"

    dropped_n <- sum(is.na(ratio_df$ratio))
    if (dropped_n > 0) {
        warning(
            dropped_n,
            " patient(s) have zero ", denominator_cluster,
            " cells; ratio set to NA for those rows."
        )
    }

    ratio_dfout <<- ratio_df

    y_breaks_fun <- function(x) {
        sort(unique(c(scales::breaks_pretty(n = 5)(x), 1)))
    }

    p <- ggplot(ratio_df, aes_string(x = group_by, y = "ratio")) +
        labs(x = NULL, y = paste0(numerator_cluster, "/", denominator_cluster, " ratio")) +
        theme_bw() +
        theme(
            panel.grid        = element_blank(),
            axis.text         = element_text(color = "black"),
            aspect.ratio      = facet_ratio,
            axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.ticks        = element_line(color = "black"),
            panel.border      = element_rect(color = "black", fill = NA, size = 0.5),
            panel.spacing     = unit(panel_spacing, "lines"),
            legend.key.height = unit(1.5, "lines")
        ) +
        # Unity reference line (ratio == 1). Drawn before boxes so it sits beneath.
        geom_hline(yintercept = 1, linetype = "dotted", color = "grey40", linewidth = 0.5) +
        geom_boxplot(
            aes_string(fill = group_by),
            color = "grey16", position = position_dodge(),
            size = 0.5, alpha = 0.8, outlier.color = NA, show.legend = TRUE
        )

    if (!is.null(shape_by)) {
        p <- p + geom_quasirandom(
            aes_string(shape = shape_by),
            size = point_size, width = 0.2
        )
    } else {
        p <- p + geom_quasirandom(
            fill = "grey84", size = point_size, width = 0.2,
            shape = 21, alpha = 0.8
        )
    }

    if (log) {
        p <- p +
            scale_y_continuous(
                trans = "log10", limits = c(miny, maxy),
                breaks = c(0.01, 0.1, 1, 10, 100),
                labels = c(0.01, 0.1, 1, 10, 100)
            ) +
            annotation_logticks(base = 10, sides = "l", outside = TRUE) +
            coord_cartesian(clip = "off") +
            scale_size_area(max_size = 15) +
            theme(axis.text.y = element_text(margin = margin(r = 8)))
    } else {
        p <- p +
            scale_y_continuous(limits = c(0, NA), breaks = y_breaks_fun) +
            expand_limits(y = 1) +
            coord_cartesian(clip = "off") +
            scale_size_area(max_size = 15)
    }

    if (show_stats) {
        stat_df <- ratio_df %>% dplyr::filter(!is.na(ratio))

        # --- NPC (Prism-style) path: single-panel via dummy facet col ---
        if (identical(stats_style, "npc")) {
            # add_pvalue_npc expects a facet variable. Inject a constant dummy.
            stat_df$.dummy_facet <- "1"

            if (is.null(external_stats)) {
                npc_stat <- tryCatch(
                    rstatix::dunn_test(stat_df,
                                       as.formula(paste("ratio ~", group_by)),
                                       p.adjust.method = "holm"),
                    error = function(e) NULL
                )
                if (!is.null(npc_stat) && nrow(npc_stat) > 0) {
                    npc_stat$.dummy_facet <- "1"
                }
            } else {
                npc_stat <- external_stats
                npc_stat$.dummy_facet <- "1"
            }

            grp_levels <- if (is.factor(stat_df[[group_by]])) {
                levels(droplevels(stat_df[[group_by]]))
            } else {
                sort(unique(as.character(stat_df[[group_by]])))
            }

            panel_max_df <- data.frame(
                .dummy_facet = "1",
                max_val = max(stat_df$ratio, na.rm = TRUE),
                stringsAsFactors = FALSE
            )
            if (panel_max_df$max_val <= 0 || is.na(panel_max_df$max_val)) {
                panel_max_df$max_val <- 1
            }

            statout <<- npc_stat
            if (!is.null(npc_stat) && nrow(npc_stat) > 0) {
                p <- p + add_pvalue_npc(
                    stat.test              = npc_stat,
                    panel_max              = panel_max_df,
                    facet_var              = ".dummy_facet",
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
            dummy_stat <- stat_df %>%
                rstatix::wilcox_test(as.formula(paste("ratio ~", group_by)), paired = FALSE) %>%
                rstatix::add_y_position(step.increase = step_increase)

            external_stats <- external_stats %>%
                dplyr::left_join(
                    dummy_stat %>% dplyr::select(group1, group2, y.position),
                    by = c("group1", "group2")
                )
            stat.test <- external_stats
        } else {
            stat.test <- stat_df %>%
                rstatix::dunn_test(as.formula(paste("ratio ~", group_by)), p.adjust.method = "holm") %>%
                rstatix::add_y_position(step.increase = step_increase)
        }

        statout <<- stat.test

        p <- p + ggpubr::stat_pvalue_manual(
            stat.test, label = "p.adj.signif",
            tip.length = 0.01, hide.ns = hide_ns, size = 5
        )
    }

    return(p)
}
