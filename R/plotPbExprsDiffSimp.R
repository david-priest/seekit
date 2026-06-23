# plotPbExprsDiffSimp — migrated from CyTOF nBass_helpers.R into seekit (CATALYST-free).
# 2026-06-10: lifted verbatim, de-CATALYST'd (.wl_* internals, namespace hack removed).

plotPbExprsDiffSimp <- function (
    x,
    scales = "free_y",
    conditions = NULL,
    color_by = NULL,
    excluded_clusters = NULL,
    axes = "all",
    features = "state",
    assay = "exprs",
    fun = c("median", "mean", "sum"),
    point_size = 1,
    clusters_order = NULL,
    textsize = 14,
    panel_spacing = 2,
    show_stats = TRUE,
    line = FALSE,
    hide_ns = TRUE,
    facet_by = NULL,
    x_group = "condition",
    mean_or_med = "median",
    group_by = x_group,
    shape_by = NULL,
    size_by = FALSE,
    point_alpha = 1,
    stat_size = 4,
    facet_ratio = 1.5,
    geom = c("boxes","bar"),
    jitter = TRUE,
    geom_type = c("quasirandom", "point"),
    group1keep = NULL,
    nudge = 0,
    step_increase = 0.1,
    hide_x_labels = FALSE,
    label_parse = FALSE,
    merging_col = NULL,
    swap = FALSE,
    external_stats = NULL,
    average_samples = FALSE
) {
    fun <- match.arg(fun)
    geom <- match.arg(geom)
    geom_type <- match.arg(geom_type)
    stopifnot(is.logical(jitter), length(jitter) == 1)

    .wl_check_sce(x)
    .wl_check_assay(x, assay)
    .wl_check_cd_factor(x, x_group)
    .wl_check_cd_factor(x, group_by)

    x <- x[.wl_get_features(x, features), ]

    by <- "sample_id"
    ms <- .wl_agg(x, by, fun, assay)
    df <- melt(ms, varnames = c("antigen", by))

    DFCOND <<- df

    i <- match(df$sample_id, x$sample_id)
    j <- setdiff(names(colData(x)), c(names(df)))
    df <- cbind(df, colData(x)[i, j, drop = FALSE])

    if (size_by) {
        size_by <- "n_cells"
        df$n_cells <- table(.wl_sample_ids(x))[df$sample_id]
    } else {
        size_by <- NULL
    }

    if (!is.null(conditions)) {
        df <- df[df[[x_group]] %in% conditions, ]
        df[[x_group]] <- factor(df[[x_group]], levels = conditions, ordered = TRUE)
    }

    if (!is.null(excluded_clusters)) {
        df <- df[!df$cluster_id %in% excluded_clusters, ]
    }

    dfout <<- df

    if (average_samples) {
        df <- df %>%
            group_by(across(all_of(c(color_by, x_group, "antigen")))) %>%
            summarise(value = mean(value, na.rm = TRUE), .groups = 'drop')
    }

    if (nrow(df) == 0) {
        stop("No data available after filtering. Please adjust the parameters.")
    }

    dfoutout <<- df

    if (show_stats) {
        if (!is.null(external_stats)) {
            external_stats <- external_stats %>% filter(antigen %in% features)
            dummy_stat <- df %>% group_by(antigen, !!sym(facet_by)) %>% wilcox_test(as.formula(paste("value ~", x_group)))
            dummy_stat <- dummy_stat %>% add_y_position(scales = "free", step.increase = step_increase)
            external_stats <- external_stats %>%
                left_join(dummy_stat %>% select(antigen, !!sym(facet_by), group1, group2, y.position),
                          by = c("antigen", facet_by, "group1", "group2"))
            stat.test <- external_stats
        } else {
            stat.test <- df %>%
                group_by(antigen, !!sym(facet_by)) %>%
                wilcox_test(as.formula(paste("value ~", x_group))) %>%
                add_y_position(scales = "free", step.increase = step_increase)
        }

        if (!is.null(group1keep)) {
            stat.test <- stat.test[stat.test$group1 %in% group1keep,]
            stat.test <- stat.test %>% add_y_position(scales = "free", step.increase = step_increase)
        }

        pbCondStats <<- stat.test
    }

    if (mean_or_med == "median") {
        fun = "median"
    } else {
        fun = "mean"
    }

    p <- ggplot(df, aes_string(x = group_by, y = "value", fill = color_by)) +
        labs(y = paste(fun, ifelse(assay == "exprs", "expression", assay))) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            text = element_text(size = textsize),
            strip.text = element_text(size = textsize),
            strip.background = element_rect(fill = NA, color = NA),
            strip.placement = "outside",
            legend.text = element_text(size = textsize),
            legend.title = element_text(size = textsize),
            aspect.ratio = facet_ratio,
            axis.ticks = element_line(color = "black"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.spacing = unit(panel_spacing, "lines"),
            axis.text = element_text(color = "black", size = textsize),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title = element_text(size = textsize)
        )

    if (geom == "boxes") {
      p <- p + geom_boxplot(aes_string(group = x_group), fill = "grey84", color = "black", width = 0.75, linewidth = 0.3, alpha = 0.9, outlier.color = NA, show.legend = FALSE)
    } else if (geom == "bar") {
      p <- p + stat_summary(fun = fun, geom = "bar", aes_string(group = x_group), fill = "grey84", color = "black", position = position_dodge(), alpha = 1)
    }

    # Lines linking each condition's points across the day axis.
    if (line == TRUE) {
      p <- p + geom_line(aes_string(group = color_by), linewidth = 0.4, alpha = 0.6)
    }

    if (geom_type == "quasirandom") {
      if (!is.null(shape_by)) {
        p <- p + geom_quasirandom(aes_string(group = x_group, color = color_by, shape = shape_by), width = 0.4, size = point_size, stroke = 0.5, show.legend = TRUE, alpha = point_alpha)
      } else {
        p <- p + geom_quasirandom(aes_string(group = x_group, fill = color_by), width = 0.4, size = point_size, shape = 21, color = "black", stroke = 0.5, show.legend = TRUE, alpha = point_alpha)
      }
    } else if (geom_type == "point") {
      if (!is.null(shape_by)) {
        p <- p + geom_point(aes_string(group = x_group, color = color_by, shape = shape_by), size = point_size, stroke = 0.5, show.legend = TRUE, alpha = point_alpha)
      } else {
        p <- p + geom_point(aes_string(group = x_group, fill = color_by), size = point_size, shape = 21, color = "black", stroke = 0.5, show.legend = TRUE, alpha = point_alpha)
      }
    }

    if (hide_x_labels) {
      p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    }

    if (show_stats) {
      p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.025, hide.ns = hide_ns, size = stat_size, bracket.nudge.y = nudge)
    }

    facet_formula <- if (is.null(facet_by)) {
        as.formula("~ antigen")
    } else {
        as.formula(paste(facet_by, "~ antigen"))
    }

    p <- p + facet_grid2(facet_formula, scales = scales, axes = axes, independent = "y")

    return(p)
}
