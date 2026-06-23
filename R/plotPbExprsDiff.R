library(SingleCellExperiment)
library(ggplot2)
library(ggbeeswarm)
library(scales)
library(dplyr)
library(rstatix)
library(rlang)
library(ggh4x)

plotPbExprsDiff <- function (
    x,
    k = "meta20",
    scales = "free_y",
    conditions = NULL,
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
    hide_ns = TRUE,
    facet_by = c("antigen", "cluster_id"),
    color_by = "condition",
    mean_or_med = "median",
    group_by = color_by,
    shape_by = NULL,
    size_by = FALSE,
    stat_size = 4,
    facet_ratio = 1.5,
    geom = c("boxes","bar"),
    jitter = TRUE,
    group1keep = NULL,
    nudge = 0,
    step_increase = 0,
    hide_x_labels = FALSE,
    label_parse = FALSE,
    merging_col = NULL,
    swap = FALSE,
    external_stats = NULL  # New parameter for external stats dataframe
) {
    fun <- match.arg(fun)
    geom <- match.arg(geom)
    facet_by <- match.arg(facet_by)
    stopifnot(is.logical(jitter), length(jitter) == 1)

    # Use the merging column from colData if provided, or if k is a colData column
    if (!is.null(merging_col) || k %in% names(colData(x))) {
        cluster_ids <- factor(x[[k]])
    } else {
        # If no merging_col, proceed with normal clustering logic
        .wl_check_sce(x)
        k <- .wl_check_k(x, k)
        cluster_ids <- .wl_cluster_ids(x, k)
    }

    # Allow plotting from column data
    .wl_check_assay(x, assay)
    .wl_check_cd_factor(x, color_by)
    .wl_check_cd_factor(x, group_by)

    # Retrieve the features
    x <- x[.wl_get_features(x, features), ]
    x$cluster_id <- cluster_ids
    by <- c("cluster_id", "sample_id")
    ms <- .wl_agg(x, by, fun, assay)
    df <- melt(ms, varnames = c("antigen", by[length(by)]))

    DFCOND <<- df

    if (!is.null(cluster_ids)) {
        df$cluster_id <- df$L1
    }

    i <- match(df$sample_id, x$sample_id)
    j <- setdiff(names(colData(x)), c(names(df), "cluster_id"))
    df <- cbind(df, colData(x)[i, j, drop = FALSE])

    ncs <- table(as.list(colData(x)[by]))
    ncs <- rep(c(t(ncs)), each = nrow(x))

    if (size_by) {
        size_by <- "n_cells"
        df$n_cells <- ncs
    } else {
        size_by <- NULL
    }

    df <- df[ncs > 0, , drop = FALSE]

    if (!is.null(conditions)) {
        df <- df[df[[color_by]] %in% conditions, ]
        df[[color_by]] <- factor(df[[color_by]], levels = conditions, ordered = TRUE)
    }

    # Apply excluded clusters if needed
    if (!is.null(excluded_clusters)) {
        df <- df[!df$cluster_id %in% excluded_clusters, ]
    }

    if (show_stats) {
        if (!is.null(external_stats)) {
            # Filter external stats to include only the antigens that are plotted
            external_stats <- external_stats %>% filter(antigen %in% features)

            # Use external stats if provided
            dummy_stat <- df %>% group_by(antigen, cluster_id) %>% wilcox_test(as.formula(paste("value ~", color_by)))
            dummy_stat <- dummy_stat %>% add_y_position(scales = "free", step.increase = step_increase)

            # Merge y.position from dummy_stat to external_stats
            external_stats <- external_stats %>%
                left_join(dummy_stat %>% select(antigen, cluster_id, group1, group2, y.position),
                          by = c("antigen", "cluster_id", "group1", "group2"))

            stat.test <- external_stats
        } else {
            # Calculate stats internally using wilcox_test (change to Dunn?)
            stat.test <- df %>%
                group_by(antigen, cluster_id) %>%
                wilcox_test(as.formula(paste("value ~", color_by))) %>%
                add_y_position(scales = "free", step.increase = step_increase)
        }

        if (!is.null(group1keep)) {
            stat.test <- stat.test[stat.test$group1 %in% group1keep,]
            stat.test <- stat.test %>% add_y_position(scales = "free", step.increase = step_increase)
        }

        # Save unfiltered stats globally
        pbCondStats <<- stat.test
    }

    if (mean_or_med == "median") {
        fun = "median"
    } else {
        fun = "mean"
    }

    # Initialize the ggplot
    p <- ggplot(df, aes_string(x = group_by, y = "value")) +
        labs(y = paste(fun, ifelse(assay == "exprs", "expression", assay))) + # Plots "median expression"
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            text = element_text(size = textsize),
            strip.text = element_text(size = textsize),
            strip.background = element_rect(fill = NA, color = NA),
            strip.placement = "outside",  # Move facet labels outside
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

    # Add geom based on selection
    if (geom == "boxes") {
        p <- p + geom_boxplot(aes_string(fill = color_by), color = "black", width = 0.75, linewidth = 0.3, alpha = 0.9, outlier.color = NA, show.legend = T)
    } else if (geom == "bar") {
        p <- p + stat_summary(fun = fun, geom = "bar", aes_string(fill = color_by), color = "black", position = position_dodge(), alpha = 1)
    }

    # Define the position adjustment for jitter
    position_adjustment <- position_quasirandom(width = 0.2)

    p <- p + geom_quasirandom(width = 0.2, size = point_size, shape = 21, fill = "grey84", color = "black", stroke = 0.5, show.legend = TRUE, alpha = 0.8)

    if (hide_x_labels) {
        p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    }

    # Add stats
    if (show_stats) {
        p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.025, hide.ns = hide_ns, size = stat_size, bracket.nudge.y = nudge)
    }

    # Label parsing option for superscripting
    if (label_parse) {
        if (!swap) {
            p <- p + facet_grid2(cluster_id ~ antigen, scales = scales, axes = axes, labeller = label_parsed)
        } else {
            p <- p + facet_grid2(antigen ~ cluster_id, scales = scales, axes = axes, labeller = label_parsed)
        }
    } else {
        if (!swap) {
            p <- p + facet_grid2(cluster_id ~ antigen, scales = scales, axes = axes, independent = "y") #, independent = "y" was used before
        } else {
            p <- p + facet_grid2(antigen ~ cluster_id, scales = scales, axes = axes, independent = "y")
        }
    }

    return(p)
}

