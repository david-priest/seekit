# plotAbundancesDiff.R — CATALYST-free rewrite promoted from dev/catalyst_quarantine/.
# Body verbatim from the project-vendored copy; only the CATALYST namespace
# shims (CATALYST:::.* internals, bare accessors, the asNamespace('CATALYST')
# hack) were rewritten to the package's .wl_* internals (R/wl_internals.R).
# 2026-06 seekit migration of the CMV CyTOF pipeline.

plotAbundancesDiff <- function(
    x,
    k             = "meta20",
    by            = c("sample_id", "cluster_id"),
    group_by      = "condition",
    shape_by      = NULL,
    col_clust     = TRUE,
    n_cols        = 4,
    log           = FALSE,
    miny          = 0.01,
    maxy          = NA,
    panel_spacing = 1,
    point_size    = 2,
    facet_ratio   = 1.5,
    label_parse   = FALSE,
    distance      = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
    linkage       = c("average", "ward.D", "single", "complete", "mcquitty", "median", "centroid", "ward.D2"),
    k_pal         = .wl_cluster_cols,
    clusters_order = NULL,
    merging_col   = NULL,   # kept for backwards compatibility; auto-detection now used
    step_increase = 0,
    external_stats = NULL,
    show_stats    = FALSE,
    hide_ns       = FALSE,         # TRUE hides non-significant brackets (Cell Press friendly)
    stats_style   = c("data", "npc"),   # "data" = legacy stat_pvalue_manual in data space (per-facet inconsistent)
                                        # "npc"  = brackets in the gutter above each panel via add_pvalue_npc (uniform spacing across facets)
    npc_group_levels = NULL,       # only used when stats_style = "npc"; factor levels of group_by in plot order
    npc_bracket_top_offset = 0.20, # top bracket sits at panel_max * (1 + this); higher = more clearance from data
    npc_bracket_step       = 0.10, # successive brackets step by panel_max * this
    npc_label_size         = 4,
    npc_richtext           = FALSE,
    npc_asterisk_pt_multiplier = 1.4,
    npc_panel_spacing_y_lines  = 3,    # vertical gap between facet rows (room for brackets above panels)
    npc_plot_margin_top_pt     = 30,
    npc_nature_style       = FALSE     # TRUE: drop panel borders + use L-shape axes (Nature/Immunity look)
) {
    stats_style <- match.arg(stats_style)
    by      <- match.arg(by)
    linkage <- match.arg(linkage)
    distance <- match.arg(distance)

    .wl_check_sce(x, TRUE)

    # Auto-detect whether k is a colData column or a cluster_codes key.
    # merging_col is kept for backwards compatibility but is no longer needed.
    if (k %in% names(colData(x))) {
        cluster_ids <- x[[k]]
    } else {
        k <- .wl_check_k(x, k)
        cluster_ids <- .wl_cluster_ids(x, k)
    }

    .wl_check_cd_factor(x, group_by)
    .wl_check_cd_factor(x, shape_by)
    .wl_check_pal(k_pal)

    stopifnot(is.logical(col_clust), length(col_clust) == 1)

    shapes <- .wl_get_shapes(x, shape_by)
    if (is.null(shapes)) shape_by <- NULL

    if (by == "sample_id") {
        nk <- nlevels(factor(cluster_ids))
        if (length(k_pal) < nk) k_pal <- colorRampPalette(k_pal)(nk)
    }

    ns  <- table(cluster_id = cluster_ids, sample_id = .wl_sample_ids(x))
    fq  <- prop.table(ns, 2) * 100
    df  <- as.data.frame(fq)
    m   <- match(df$sample_id, x$sample_id)

    for (i in c(shape_by, group_by)) df[[i]] <- x[[i]][m]

    if (by == "sample_id" && col_clust && length(unique(df$sample_id)) > 1) {
        d  <- dist(t(fq), distance)
        h  <- hclust(d, linkage)
        o  <- colnames(fq)[h$order]
        df$sample_id <- factor(df$sample_id, o)
    }

    if (!is.null(clusters_order)) {
        df$cluster_id <- factor(df$cluster_id, levels = clusters_order)
    }

    if (log == TRUE) {
        df$Freq <- df$Freq + 0.01
    }

    dfout <<- df

    maxy_values <- df %>%
        dplyr::group_by(cluster_id) %>%
        dplyr::summarise(maxy = ceiling(max(Freq, na.rm = TRUE) * 1.1), .groups = "drop")

    p <- ggplot(df, aes_string(y = "Freq")) +
        labs(x = NULL, y = "Proportion [%]") +
        theme_bw() +
        theme(
            panel.grid       = element_blank(),
            strip.text       = element_text(face = "bold"),
            strip.background = element_rect(fill = NA, color = NA),
            axis.text        = element_text(color = "black"),
            aspect.ratio     = facet_ratio,
            axis.text.x      = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.ticks       = element_line(color = "black"),
            panel.border     = element_rect(color = "black", fill = NA, size = 0.5),
            panel.spacing    = unit(panel_spacing, "lines"),
            legend.key.height = unit(1.5, "lines")
        )

    if (label_parse) {
        p <- p +
            facet_wrap(~cluster_id, scales = "free", ncol = n_cols, labeller = label_parsed) +
            geom_boxplot(
                aes_string(x = group_by, fill = group_by),
                color = "grey16", position = position_dodge(),
                size = 0.5, alpha = 0.8, outlier.color = NA, show.legend = TRUE
            )
    } else {
        p <- p +
            facet_wrap(~cluster_id, scales = "free", ncol = n_cols,
                       labeller = labeller(cluster_id = label_value)) +
            geom_boxplot(
                aes_string(x = group_by, fill = group_by),
                color = "grey16", position = position_dodge(),
                size = 0.5, alpha = 0.8, outlier.color = NA, show.legend = TRUE
            )
    }

    if (!is.null(shape_by)) {
        p <- p + geom_quasirandom(
            aes_string(x = group_by, shape = shape_by),
            size = point_size, width = 0.2
        )
    } else {
        p <- p + geom_quasirandom(
            aes_string(x = group_by),
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
        # When npc_nature_style = TRUE, use ggh4x::facet_wrap2(axes = "all_x")
        # so x-axes (ticks + labels) appear on every panel, not just bottom row.
        # This matches Nature/Immunity style where multi-row figures keep axes
        # on each row.
        facet_layer <- if (isTRUE(npc_nature_style) && requireNamespace("ggh4x", quietly = TRUE)) {
            # ggh4x::facet_wrap2 axes arg: "margins" (default), "x", "y", "all"
            ggh4x::facet_wrap2(~cluster_id, scales = "free_y", ncol = n_cols,
                               axes = "x",
                               labeller = labeller(cluster_id = label_value))
        } else {
            facet_wrap(~cluster_id, scales = "free_y", ncol = n_cols,
                       labeller = labeller(cluster_id = label_value))
        }
        p <- p +
            scale_y_continuous(limits = c(0, NA)) +
            coord_cartesian(clip = "off") +
            scale_size_area(max_size = 15) +
            facet_layer
    }

    p <- p + geom_blank(data = maxy_values, aes(y = maxy))

    if (show_stats) {
        # --- NPC stats path: panel-relative brackets, uniform across facets ---
        # Bypasses stat_pvalue_manual entirely; add_pvalue_npc handles hide_ns,
        # ranking, axis-extension, and bracket rendering itself.
        if (identical(stats_style, "npc")) {
            if (is.null(external_stats)) {
                external_stats <- df %>%
                    dplyr::group_by(cluster_id) %>%
                    rstatix::dunn_test(as.formula(paste("Freq ~", group_by)), p.adjust.method = "holm") %>%
                    dplyr::ungroup()
            }
            panel_max_df <- df %>%
                dplyr::group_by(cluster_id) %>%
                dplyr::summarise(max_val = max(Freq, na.rm = TRUE), .groups = "drop") %>%
                dplyr::mutate(max_val = ifelse(max_val == 0 | is.na(max_val), 1, max_val))

            grp_levels <- if (!is.null(npc_group_levels)) {
                npc_group_levels
            } else if (is.factor(df[[group_by]])) {
                levels(droplevels(df[[group_by]]))
            } else {
                sort(unique(as.character(df[[group_by]])))
            }

            statout <<- external_stats

            p <- p + add_pvalue_npc(
                stat.test              = external_stats,
                panel_max              = panel_max_df,
                facet_var              = "cluster_id",
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

            return(p)
        }

        # --- Legacy data-space stats path (default) -----------------------
        if (!is.null(external_stats)) {
            dummy_stat <- df %>%
                dplyr::group_by(cluster_id) %>%
                rstatix::wilcox_test(as.formula(paste("Freq ~", group_by)), paired = FALSE) %>%
                rstatix::add_y_position(scales = "free", step.increase = step_increase)

            external_stats <- external_stats %>%
                dplyr::left_join(
                    dummy_stat %>% dplyr::select(cluster_id, group1, group2, y.position),
                    by = c("cluster_id", "group1", "group2")
                )
            stat.test <- external_stats
        } else {
            stat.test <- df %>%
                dplyr::group_by(cluster_id) %>%
                rstatix::dunn_test(as.formula(paste("Freq ~", group_by)), p.adjust.method = "holm") %>%
                rstatix::add_y_position(scales = "free", step.increase = step_increase)
        }

        # Compact significant-only brackets when hide_ns = TRUE.
        # ggpubr::stat_pvalue_manual(hide.ns = TRUE) just skips drawing ns rows
        # but leaves their y.position slots empty -> visible gaps. We instead
        # drop ns rows and re-rank y.position so brackets stack tight.
        if (isTRUE(hide_ns) && "p.adj.signif" %in% colnames(stat.test)) {
            stat.test <- stat.test %>% dplyr::filter(p.adj.signif != "ns")

            if (nrow(stat.test) > 0) {
                max_y_df <- df %>%
                    dplyr::group_by(cluster_id) %>%
                    dplyr::summarise(max_val = max(Freq, na.rm = TRUE), .groups = "drop") %>%
                    dplyr::mutate(max_val = ifelse(max_val == 0 | is.na(max_val), 1, max_val))

                eff_step <- if (step_increase == 0) 0.1 else step_increase

                stat.test <- stat.test %>%
                    dplyr::select(-dplyr::any_of("y.position")) %>%
                    dplyr::left_join(max_y_df, by = "cluster_id") %>%
                    dplyr::group_by(cluster_id) %>%
                    dplyr::mutate(
                        y.position = max_val * (1 + eff_step * dplyr::row_number())
                    ) %>%
                    dplyr::ungroup() %>%
                    dplyr::select(-max_val)
            }
        }

        statout <<- stat.test

        if (!is.null(stat.test) && nrow(stat.test) > 0) {
            p <- p + ggpubr::stat_pvalue_manual(
                stat.test, label = "p.adj.signif",
                tip.length = 0.01,
                # rows already filtered above when hide_ns = TRUE, so hide.ns=FALSE here
                hide.ns = FALSE, size = 5
            )
        }
    }

    return(p)
}

