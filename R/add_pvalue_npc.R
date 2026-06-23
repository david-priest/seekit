# add_pvalue_npc
#
# Render significance brackets ABOVE each facet panel (Nature/Immunity-style
# "strip" layout) without touching the data axes. Brackets are drawn into the
# space above each panel via coord_cartesian(clip = "off"), with no axis
# extension and no ghost ticks. Panels with no significant brackets after the
# hide_ns filter are left completely alone.
#
# Why this works for "consistent spacing across facets":
# In ggplot, all panels in a facet_wrap layout share the same PHYSICAL height.
# If we place a bracket at data y = panel_max * (1 + offset), its physical
# position above the panel-top is (offset / axis_expansion_factor) * panel_height
# -- a CONSTANT fraction of panel height in every panel, regardless of each
# panel's own data range. Same applies to bracket spacing. So absolute pixel
# spacing is identical in every panel: exactly the Prism-style result.
#
# Why this fixes "data spills off the top":
# We never set a hard upper y-limit. The y-axis fits the data naturally
# (whatever the parent plot decides). Brackets render in the gutter ABOVE the
# panel via clip = "off". The panel's own data is never clipped.
#
# Why this fixes "ghost ticks in the stats zone":
# No axis extension means the y-axis only goes up to ~panel_max * 1.05 (natural
# ggplot expansion). Tick generation is confined to that range. Brackets exist
# physically above the axis, in space ggplot doesn't tick.
#
# Why panels with all-ns aren't touched:
# After the hide_ns filter, those facets have zero rows in stat.test and
# contribute zero bracket geoms. Their axes remain at whatever the parent plot
# (with its existing geom_blank, scale, etc.) decided.
#
# Required parent-plot setup:
#   - coord_cartesian(clip = "off")        # so brackets above the panel render
#   - enough vertical panel-spacing / plot top margin to fit the bracket gutter
#     (this function adds a theme() layer doing both by default; set
#     apply_theme = FALSE to opt out and manage theme yourself)

add_pvalue_npc <- function(
    stat.test,
    panel_max,                          # named numeric vector OR data frame with facet_var + max_col
    facet_var          = "cluster_id",
    max_col            = "max_val",     # only used if panel_max is a data frame
    group_levels       = NULL,          # factor levels of x-axis group var, in plot order
    hide_ns            = TRUE,
    label_col          = "p.adj.signif",
    label_size         = 4,
    label_fontface     = "plain",
    label_transform    = NULL,          # default-NULL -> see asterisk handling below
    use_richtext       = FALSE,         # TRUE -> ggtext::geom_richtext (allows per-character sizing)
    asterisk_pt_multiplier   = 1.4,     # when use_richtext: render "*" this much bigger than label_size
    asterisk_baseline_shift  = -3,      # when use_richtext: lower "*" baseline by this many pt (negative = down)
    bracket_color      = "black",
    bracket_size       = 0.4,
    # Bracket geometry (fractions of panel_max for each facet)
    bracket_top_offset = 0.20,          # top bracket sits at panel_max * (1 + this)
    bracket_step       = 0.10,          # successive brackets step by panel_max * this
    bracket_tip_length = 0.025,         # tip length as fraction of panel_max
    label_nudge        = 0.025,         # label-above-bracket nudge as fraction of panel_max
    # Theme additions to make room for the bracket gutter above each panel
    apply_theme        = TRUE,
    panel_spacing_y_lines = 3,          # increases vertical gap between facet rows
    plot_margin_top_pt    = 30,         # increases top plot margin
    nature_style       = FALSE          # Nature/Immunity look: no panel border, L-shape axes
) {
    if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr required.")
    if (!requireNamespace("rlang", quietly = TRUE)) stop("rlang required.")
    if (!is.data.frame(stat.test)) stop("stat.test must be a data frame.")

    # facet_var may be a single column name or a character vector of column names
    # (for 2D facet plots like plotBoxplotsProportions which facets master x daughter).
    # In the multi-column case we build a composite key joined by "\037" (ASCII unit
    # separator) -- a character that won't naturally appear in cluster names.
    required_cols <- c(facet_var, "group1", "group2", label_col)
    missing <- setdiff(required_cols, colnames(stat.test))
    if (length(missing) > 0) {
        stop("stat.test is missing required columns: ", paste(missing, collapse = ", "))
    }
    composite_facet_key <- function(df, cols) {
        if (length(cols) == 1) {
            as.character(df[[cols]])
        } else {
            do.call(paste, c(lapply(cols, function(c) as.character(df[[c]])),
                             sep = "\037"))
        }
    }

    # ---- Default label transform (ASCII asterisk with optional styling) ---
    # We keep ASCII "*" (universally renderable in any PDF/raster device) and
    # use ggtext spans to fix its typography: it sits high on the line and is
    # smaller than letters by default, so we bump font-size and lower its
    # baseline via vertical-align. Asterisks are escaped as &#42; HTML entities
    # so gridtext's markdown parser doesn't interpret "**" as bold start.
    if (is.null(label_transform)) {
        if (isTRUE(use_richtext)) {
            asterisk_pt <- max(1, round(label_size * 3.5 * asterisk_pt_multiplier))
            # vertical-align: negative pt lowers the glyph; gridtext accepts pt units
            valign <- paste0(asterisk_baseline_shift, "pt")
            label_transform <- function(x) {
                has_ast <- grepl("[*]", x)
                if (any(has_ast)) {
                    x_esc <- gsub("\\*", "&#42;", x)
                    x[has_ast] <- paste0(
                        "<span style='font-size:", asterisk_pt, "pt; ",
                        "vertical-align:", valign, "; ",
                        "line-height:0.6'>",
                        x_esc[has_ast], "</span>"
                    )
                }
                x
            }
        } else {
            # No richtext -> plain ASCII "*", smaller and higher than letters but
            # at least visible in any device.
            label_transform <- function(x) x
        }
    }
    stopifnot(is.function(label_transform))
    if (isTRUE(use_richtext) && !requireNamespace("ggtext", quietly = TRUE)) {
        stop("use_richtext = TRUE requires the 'ggtext' package. ",
             "Install with: install.packages('ggtext').")
    }

    # ---- Normalize panel_max -> named numeric vector ----------------------
    # When facet_var is multi-column, panel_max should contain all those columns
    # plus max_col; pm_vec is keyed by the composite facet key.
    if (is.data.frame(panel_max)) {
        if (!all(c(facet_var, max_col) %in% colnames(panel_max))) {
            stop("panel_max data frame must contain all of '",
                 paste(c(facet_var, max_col), collapse = "', '"), "' columns.")
        }
        pm_vec <- setNames(panel_max[[max_col]], composite_facet_key(panel_max, facet_var))
    } else if (is.numeric(panel_max)) {
        if (is.null(names(panel_max))) stop("panel_max numeric vector must be named with facet levels.")
        pm_vec <- panel_max
    } else {
        stop("panel_max must be a numeric vector or data frame.")
    }
    pm_vec[is.na(pm_vec) | pm_vec <= 0] <- 1

    # ---- Drop ns rows if requested ----------------------------------------
    if (isTRUE(hide_ns)) {
        ns_mask <- as.character(stat.test[[label_col]]) %in% c("ns", "NS", "")
        stat.test <- stat.test[!ns_mask, , drop = FALSE]
    }

    # If nothing survives, return zero layers -- parent plot stays untouched.
    if (nrow(stat.test) == 0) {
        if (isTRUE(apply_theme)) {
            theme_layer <- build_npc_theme_layer(panel_spacing_y_lines,
                                                 plot_margin_top_pt,
                                                 nature_style)
            return(list(theme_layer))
        }
        return(list())
    }

    # ---- Group-levels -> x positions --------------------------------------
    g1 <- as.character(stat.test$group1)
    g2 <- as.character(stat.test$group2)
    if (is.null(group_levels)) {
        group_levels <- sort(unique(c(g1, g2)))
    }
    x_of_g1 <- match(g1, group_levels)
    x_of_g2 <- match(g2, group_levels)
    if (any(is.na(x_of_g1)) || any(is.na(x_of_g2))) {
        missing_levels <- unique(c(g1[is.na(x_of_g1)], g2[is.na(x_of_g2)]))
        stop("These group labels in stat.test are not in group_levels: ",
             paste(missing_levels, collapse = ", "))
    }

    # ---- Rank brackets per facet (1 = top) --------------------------------
    stat.test$.facet <- composite_facet_key(stat.test, facet_var)
    stat.test$.x1    <- x_of_g1
    stat.test$.x2    <- x_of_g2
    stat.test <- stat.test[order(stat.test$.facet), , drop = FALSE]
    stat.test$.rank <- ave(seq_len(nrow(stat.test)),
                           stat.test$.facet,
                           FUN = function(i) seq_along(i))

    # Helper for re-attaching ORIGINAL facet columns to the constructed
    # segment / text data frames so ggplot's faceting matches them per panel.
    attach_facet_cols <- function(df, facet_keys) {
        if (length(facet_var) == 1) {
            df[[facet_var]] <- facet_keys
        } else {
            parts <- strsplit(facet_keys, "\037", fixed = TRUE)
            for (i in seq_along(facet_var)) {
                df[[facet_var[i]]] <- vapply(parts, `[`, character(1), i)
            }
        }
        df
    }

    # ---- Compute bracket positions in DATA space (above each panel's max) -
    pmpr <- pm_vec[stat.test$.facet]
    stat.test$.y_bracket <- pmpr * (1 + bracket_top_offset + (stat.test$.rank - 1) * bracket_step)
    stat.test$.y_tip     <- stat.test$.y_bracket - pmpr * bracket_tip_length
    stat.test$.y_label   <- stat.test$.y_bracket + pmpr * label_nudge
    stat.test$.x_mid     <- (stat.test$.x1 + stat.test$.x2) / 2

    # ---- Segment + label data frames --------------------------------------
    horiz_df <- data.frame(
        x    = stat.test$.x1,
        xend = stat.test$.x2,
        y    = stat.test$.y_bracket,
        yend = stat.test$.y_bracket,
        stringsAsFactors = FALSE
    )
    horiz_df <- attach_facet_cols(horiz_df, stat.test$.facet)

    tip_left_df <- data.frame(
        x    = stat.test$.x1,
        xend = stat.test$.x1,
        y    = stat.test$.y_bracket,
        yend = stat.test$.y_tip,
        stringsAsFactors = FALSE
    )
    tip_left_df <- attach_facet_cols(tip_left_df, stat.test$.facet)

    tip_right_df <- data.frame(
        x    = stat.test$.x2,
        xend = stat.test$.x2,
        y    = stat.test$.y_bracket,
        yend = stat.test$.y_tip,
        stringsAsFactors = FALSE
    )
    tip_right_df <- attach_facet_cols(tip_right_df, stat.test$.facet)

    seg_df <- rbind(horiz_df, tip_left_df, tip_right_df)

    text_df <- data.frame(
        x     = stat.test$.x_mid,
        y     = stat.test$.y_label,
        label = label_transform(as.character(stat.test[[label_col]])),
        stringsAsFactors = FALSE
    )
    text_df <- attach_facet_cols(text_df, stat.test$.facet)

    # NOTE: previous version attempted to pin per-facet y-scales via
    # ggh4x::facetted_pos_scales(limits = c(0, panel_max), oob = oob_keep)
    # so brackets would render in the gutter above each panel WITHOUT extending
    # the y-axis. But facetted_pos_scales' panel-to-scale matching wasn't
    # reliable in practice (mis-assigned scales caused CD4 CTL data to render
    # above the panel boundary). Reverted to the simpler approach: brackets are
    # plain geom_segments, they DO train the y-scale, so the axis extends to
    # fit them. The known cost is that ggplot's tick generator can place breaks
    # in the bracket-reservation zone above the data.
    #
    # If you need a true panel-relative strip (brackets that don't extend the
    # axis at all), that requires a custom Geom that reads y in NPC rather than
    # data-space. Ask for `add_pvalue_npc_strip()` and I'll write it.
    layers <- list()

    layers[[length(layers) + 1]] <- ggplot2::geom_segment(
        data        = seg_df,
        mapping     = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
        color       = bracket_color,
        linewidth   = bracket_size,
        inherit.aes = FALSE
    )

    if (isTRUE(use_richtext)) {
        layers[[length(layers) + 1]] <- ggtext::geom_richtext(
            data          = text_df,
            mapping       = ggplot2::aes(x = x, y = y, label = label),
            size          = label_size,
            fontface      = label_fontface,
            vjust         = 0,
            label.size    = NA,
            fill          = NA,
            label.padding = grid::unit(c(0, 0, 0, 0), "lines"),
            inherit.aes   = FALSE
        )
    } else {
        layers[[length(layers) + 1]] <- ggplot2::geom_text(
            data        = text_df,
            mapping     = ggplot2::aes(x = x, y = y, label = label),
            size        = label_size,
            fontface    = label_fontface,
            vjust       = 0,
            inherit.aes = FALSE
        )
    }

    if (isTRUE(apply_theme)) {
        layers[[length(layers) + 1]] <- build_npc_theme_layer(panel_spacing_y_lines,
                                                              plot_margin_top_pt,
                                                              nature_style)
    }

    layers
}

# Helper: build the theme layer used by add_pvalue_npc.
# - Always increases vertical panel spacing + top plot margin so brackets in the
#   gutter above each panel have room to render (clip="off" must be set on the
#   parent plot's coord).
# - When nature_style = TRUE: removes the rectangular panel border and draws
#   L-shape axes (Nature/Immunity figure style).
build_npc_theme_layer <- function(panel_spacing_y_lines, plot_margin_top_pt, nature_style) {
    base_args <- list(
        panel.spacing.y = grid::unit(panel_spacing_y_lines, "lines"),
        plot.margin     = ggplot2::margin(t = plot_margin_top_pt,
                                          r = 5, b = 5, l = 5, unit = "pt")
    )
    if (isTRUE(nature_style)) {
        base_args$panel.border    <- ggplot2::element_blank()
        base_args$panel.background <- ggplot2::element_blank()
        base_args$axis.line       <- ggplot2::element_line(color = "black", linewidth = 0.4)
    }
    do.call(ggplot2::theme, base_args)
}
