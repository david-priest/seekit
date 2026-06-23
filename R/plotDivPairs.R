# plotDivPairs — migrated from CyTOF nBass_helpers.R into seekit (CATALYST-free).
# 2026-06-10: lifted verbatim, de-CATALYST'd (.wl_* internals, namespace hack removed).

plotDivPairs <- function (
    x,
    div_col = "div",
    parse_div = TRUE,
    metric = c("mean_div", "median_div", "percent_proliferated"),
    scales = "free_y",
    conditions = NULL,
    condition_pairs = NULL,
    facet_by = NULL,
    facet = TRUE,
    excluded_clusters = NULL,
    axes = "all",
    fun = c("median", "mean", "sum"),
    point_size = 1,
    textsize = 14,
    panel_spacing = 2,
    show_stats = TRUE,
    hide_ns = TRUE,
    fill_by = NULL,
    x_group = "day",
    shape_by = NULL,
    stat_size = 4,
    facet_ratio = 1.5,
    geom = c("boxes", "bar", "median_line"),
    jitter = TRUE,
    group1keep = NULL,
    nudge = 0,
    step_increase = 0.1,
    hide_x_labels = FALSE,
    label_parse = FALSE,
    merging_col = NULL,
    swap = FALSE,
    external_stats = NULL
) {
    fun <- match.arg(fun)
    geom <- match.arg(geom)
    metric <- match.arg(metric)

    if (!"patient_id" %in% colnames(colData(x))) stop("Column 'patient_id' is not found in colData(x).")
    if (!div_col %in% colnames(colData(x))) stop("Specified 'div_col' does not exist.")

    div_values <- colData(x)[[div_col]]
    if (parse_div) div_values <- as.numeric(gsub("Div", "", div_values))
    df <- data.frame(div = div_values, colData(x))

    if (!is.null(conditions)) {
        df <- df[df[[x_group]] %in% conditions, ]
        df[[x_group]] <- factor(df[[x_group]], levels = conditions, ordered = TRUE)
    }
    if (!is.null(excluded_clusters)) df <- df[!df$cluster_id %in% excluded_clusters, ]

    if (!is.null(condition_pairs) && length(condition_pairs) > 0) {
      df <- purrr::imap_dfr(condition_pairs, function(pair, pname) {
        df %>% dplyr::filter(!!sym(fill_by) %in% pair) %>% dplyr::mutate(condition_pair = pname)
      })
      facet_by <- "condition_pair"
    }

    agg_fun <- if (metric == "mean_div") mean else if (metric == "median_div") median
    if (metric == "percent_proliferated") {
        df_agg <- df %>%
            dplyr::group_by(sample_id, !!sym(x_group), !!sym(facet_by), patient_id, !!sym(fill_by)) %>%
            dplyr::summarise(div = sum(div != 0) / dplyr::n() * 100, .groups = 'drop')
    } else {
        df_agg <- df %>%
            dplyr::group_by(sample_id, !!sym(x_group), !!sym(facet_by), patient_id, !!sym(fill_by)) %>%
            dplyr::summarise(div = agg_fun(div, na.rm = TRUE), .groups = 'drop')
    }
    if (nrow(df_agg) == 0) stop("No data available after filtering.")

    if (!is.null(fill_by) && !is.null(facet_by)) {
        df_agg$group_median <- interaction(df_agg[[x_group]], df_agg[[fill_by]], df_agg[[facet_by]])
    } else if (!is.null(fill_by)) {
        df_agg$group_median <- interaction(df_agg[[x_group]], df_agg[[fill_by]])
    } else if (!is.null(facet_by)) {
        df_agg$group_median <- interaction(df_agg[[x_group]], df_agg[[facet_by]])
    } else {
        df_agg$group_median <- df_agg[[x_group]]
    }
    dodge_width <- 0.8
    if (!is.null(fill_by) && !is.factor(df_agg[[fill_by]])) df_agg[[fill_by]] <- factor(df_agg[[fill_by]])

    p <- ggplot(df_agg, aes(x = .data[[x_group]], y = div)) +
      labs(y = dplyr::case_when(
        metric == "percent_proliferated" ~ "Percentage Proliferated",
        metric == "mean_div" ~ "Mean Division",
        TRUE ~ "Median Division")) +
      theme_bw() +
      theme(panel.grid = element_blank(), text = element_text(size = textsize),
        strip.text = element_text(size = textsize), strip.background = element_rect(fill = NA, color = NA),
        strip.placement = "outside", legend.text = element_text(size = textsize),
        legend.title = element_text(size = textsize), aspect.ratio = facet_ratio,
        axis.ticks = element_line(color = "black"), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        panel.spacing = unit(panel_spacing, "lines"), axis.text = element_text(color = "black", size = textsize),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title = element_text(size = textsize)) +
      scale_x_discrete(expand = c(0.05, 0.05))

    if (geom == "boxes") {
      p <- p + geom_boxplot(aes(fill = .data[[fill_by]], group = interaction(.data[[x_group]], .data[[fill_by]])),
        position = position_dodge(width = dodge_width), color = "black", width = 0.8, alpha = 0.8,
        outlier.color = NA, size = 0.5, show.legend = FALSE)
    } else if (geom == "bar") {
      p <- p + stat_summary(fun = match.fun(fun), geom = "bar",
        aes(fill = .data[[fill_by]], group = interaction(.data[[x_group]], .data[[fill_by]])),
        position = position_dodge(width = dodge_width), width = 0.8, color = "black", alpha = 0.8,
        size = 0.5, show.legend = FALSE)
    }

    if (!is.null(shape_by)) {
      p <- p + ggbeeswarm::geom_quasirandom(
        aes(fill = .data[[fill_by]], shape = .data[[shape_by]], group = interaction(.data[[x_group]], .data[[fill_by]])),
        color = "black", stroke = 0.5, dodge.width = dodge_width, width = 0.1, groupOnX = TRUE,
        size = point_size, alpha = 0.9, show.legend = TRUE) +
        scale_shape_manual(values = c(21, 22, 24, 25))
    } else {
      p <- p + geom_point(
        aes(fill = .data[[fill_by]], group = interaction(.data[[x_group]], .data[[fill_by]])),
        shape = 21, color = "black", stroke = 0.5,
        position = position_jitterdodge(jitter.width = 0.1, dodge.width = dodge_width),
        size = point_size, alpha = 0.9, show.legend = TRUE)
    }

    if (geom == "median_line") {
      p <- p + stat_summary(data = df_agg, aes(x = .data[[x_group]], y = div, group = group_median),
        fun.data = function(x) { m <- median(x, na.rm = TRUE); data.frame(y = m, ymin = m, ymax = m) },
        geom = "crossbar", width = 0.7, color = "black", size = 0.8, fatten = 1,
        position = position_dodge(width = dodge_width), show.legend = FALSE)
    }

    if (!is.null(fill_by)) p <- p + labs(fill = fill_by)
    if (hide_x_labels) p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    if (facet) p <- p + facet_grid(as.formula(paste("~", facet_by)), scales = scales)
    if (metric == "percent_proliferated") p <- p + scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))
    p <- p + guides(fill = guide_legend(override.aes = list(shape = 21, size = 3)))
    return(p)
}
