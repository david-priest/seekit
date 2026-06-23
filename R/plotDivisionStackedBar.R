# plotDivisionStackedBar — migrated from CyTOF nBass_helpers.R into seekit (CATALYST-free).
# 2026-06-10: lifted verbatim, de-CATALYST'd (.wl_* internals, namespace hack removed).

plotDivisionStackedBar <- function(
    x,
    div_col = "div",
    parse_div = TRUE,
    x_group = "day",
    facet_by = "condition",
    colors = NULL
) {
    if (!div_col %in% colnames(colData(x))) {
        stop("The specified 'div_col' (", div_col, ") does not exist in the colData of the input object.")
    }

    div_values <- colData(x)[[div_col]]
    if (parse_div) {
        div_values <- as.numeric(gsub("Div", "", div_values))
    }

    df <- data.frame(div = div_values, colData(x))

    # dplyr:: prefixes: seekit (loaded before this file) attaches a
    # package that masks the bare count()/group_by(), which otherwise throws
    # "Argument 'x' is not a vector: list".
    df_prop <- df %>%
        dplyr::group_by(sample_id, !!sym(x_group), !!sym(facet_by)) %>%
        dplyr::count(div) %>%
        dplyr::mutate(prop = n / sum(n)) %>%
        dplyr::ungroup()

    df_agg <- df_prop %>%
        dplyr::group_by(!!sym(x_group), !!sym(facet_by), div) %>%
        dplyr::summarise(prop = mean(prop, na.rm = TRUE), .groups = 'drop')

    p <- ggplot(df_agg, aes(x = !!sym(x_group), y = prop, fill = factor(div))) +
        geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
        labs(
            x = x_group,
            y = "Proportion of Cells",
            fill = "Division"
        ) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            strip.background = element_rect(fill = NA, color = NA),
            strip.text = element_text(size = 12),
            # Un-boxed panels: drop the full border, show just the x and y axis lines.
            panel.border = element_blank(),
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.ticks = element_line(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            legend.position = "right"
        ) +
        facet_wrap(as.formula(paste("~", facet_by)), scales = "free_x")

    if (is.null(colors)) {
        p <- p + scale_fill_brewer(palette = "Paired")
    } else {
        p <- p + scale_fill_manual(values = colors)
    }

    return(p)
}
