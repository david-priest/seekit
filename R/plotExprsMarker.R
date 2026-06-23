# plotExprsMarker.R — CATALYST-free rewrite (body is David's own; .check_sce -> .wl_check_sce, namespace hack removed).
plotExprsMarker <- function(x, markers, condition_col = "condition",
                            facet_by = NULL, conditions = NULL,
                            condition_palette = NULL, add_border = TRUE,
                            text_size = 12, panel_spacing = 0.4,
                            alpha_amount = 0.5, ncol = NULL,
                            x_limits = c(-0.5, NA), free_y_scale = FALSE,
                            free_x_scale = TRUE, all_x_axes = FALSE)
{
    .wl_check_sce(x, TRUE)

    # Validate markers
    missing_markers <- markers[!markers %in% rownames(x)]
    if (length(missing_markers) > 0) {
        stop(paste0("Markers not found in rownames of x: ",
                    paste(missing_markers, collapse = ", ")))
    }

    # Validate condition column
    if (!condition_col %in% names(colData(x))) {
        stop(paste0("'", condition_col, "' not found in colData."))
    }

    # Validate facet column if provided
    if (!is.null(facet_by) && !facet_by %in% names(colData(x))) {
        stop(paste0("'", facet_by, "' not found in colData."))
    }

    cd <- colData(x)
    es <- assay(x[markers, , drop = FALSE], "exprs")
    df <- data.frame(t(es), cd, check.names = FALSE)
    df <- reshape2::melt(df, id.vars = names(cd),
                         variable.name = "antigen", value.name = "expression")

    # Keep antigen as an ordered factor matching the input marker order
    df$antigen <- factor(df$antigen, levels = markers)

    # Rename condition column for consistent use
    df$condition <- factor(df[[condition_col]])

    # Subset/reorder conditions if requested
    if (!is.null(conditions)) {
        df <- df[df$condition %in% conditions, ]
        df$condition <- factor(df$condition, levels = conditions)
    }

    # Build fill/colour scales
    if (!is.null(condition_palette)) {
        scale_fill  <- scale_fill_manual(values  = condition_palette)
        scale_color <- scale_color_manual(values = condition_palette)
    } else {
        scale_fill  <- scale_fill_brewer(palette  = "Set1")
        scale_color <- scale_color_brewer(palette = "Set1")
    }

    # Choose facet scales from the four combinations of free x / free y
    facet_scales <- if (free_x_scale && free_y_scale)  "free"   else
                    if (free_x_scale && !free_y_scale) "free_x" else
                    if (!free_x_scale && free_y_scale) "free_y" else
                    "fixed"

    # In grid mode, fixed breaks and limits apply globally and prevent each
    # marker column from scaling independently, so we let ggplot choose both
    # automatically per panel instead.
    grid_mode <- !is.null(facet_by) && length(markers) > 1
    active_x_limits <- if (grid_mode) c(NA, NA) else x_limits
    active_x_breaks <- if (grid_mode) waiver() else seq(0, 10, by = 2)

    p <- ggplot(df, aes(x = expression, fill = condition, color = condition)) +
        geom_density(alpha = alpha_amount, linewidth = 0.4) +
        scale_fill +
        scale_color +
        scale_x_continuous(breaks = active_x_breaks,
                           limits = active_x_limits,
                           expand = c(0, 0)) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(x = "Expression", y = "Density",
             fill = condition_col, color = condition_col) +
        theme(
            strip.background  = element_blank(),
            strip.text        = element_text(face = "bold", color = "black", size = text_size),
            panel.spacing     = unit(panel_spacing, "lines"),
            plot.margin       = unit(c(0.5, 1, 0.5, 0.5), "lines"),
            panel.background  = element_rect(fill = "transparent", color = NA),
            panel.border      = if (add_border) element_rect(color = "black", fill = NA, linewidth = 0.8) else element_blank(),
            panel.grid        = element_blank(),
            axis.line.x       = element_line(color = "black"),
            axis.line.y       = element_line(color = "black"),
            axis.text         = element_text(color = "black"),
            axis.ticks        = element_line(color = "black"),
            legend.key        = element_rect(fill = "transparent"),
            text              = element_text(size = text_size)
        )

    # Faceting: grid (rows = facet_by, cols = markers) or wrap (markers only)
    if (!is.null(facet_by)) {
        # facet_grid only supports free y per-row, not per-panel. ggh4x::facet_grid2
        # with independent = "y" gives each panel its own y scale while still
        # locking x within each marker column.
        p <- p + ggh4x::facet_grid2(
            rows = vars(!!sym(facet_by)),
            cols = vars(antigen),
            scales = "free",
            independent = "y",
            axes = if (all_x_axes) "all" else "margins"
        )
    } else if (length(markers) > 1) {
        p <- p + facet_wrap(~antigen, scales = facet_scales, ncol = ncol)
    }

    p
}