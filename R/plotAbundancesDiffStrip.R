# plotAbundancesDiffStrip.R — CATALYST-free rewrite promoted from dev/catalyst_quarantine/.
# Body verbatim from the project-vendored copy; only the CATALYST namespace
# shims (CATALYST:::.* internals, bare accessors, the asNamespace('CATALYST')
# hack) were rewritten to the package's .wl_* internals (R/wl_internals.R).
# 2026-06 seekit migration of the CMV CyTOF pipeline.
# plotAbundancesDiffStrip
#
# Sibling of plotAbundancesDiff with a TRUE separated stats strip: each cluster
# panel is composed of two stacked sub-plots,
#
#   [ stats strip ]   <- brackets only, NO axes, height scales with N brackets
#   [   data plot ]   <- boxplots + points, axes intact, untouched by brackets
#
# Stacked via patchwork::plot_layout(heights = ...).
#
# Properties (per the design discussion):
#   - The data plot's y-axis line stops at the data top (no extension into a
#     "ghost stats reservation"). Ticks only inside the data range.
#   - Strip and data plot share the same x-axis groups; the strip just doesn't
#     show axis ticks/labels.
#   - When a cluster has zero significant brackets (after hide_ns filter), it
#     gets NO strip at all -- the data plot uses the full vertical height,
#     so empty-stats clusters aren't squashed.
#   - Strip height = N_brackets * strip_height_per_bracket. So a cluster with
#     6 brackets has a tall strip; one with 1 bracket gets a thin strip.
#     Bracket spacing within the strip is constant (1.0 in strip-local units),
#     producing identical visual spacing across all clusters.
#
# Returns a patchwork object; can be `+ theme(...)`'d like a ggplot and saved
# via f2()/f() (which now uses cairo_pdf for Unicode-safe output).

# Dependencies (loaded by SETUP 1 of the calling .Rmd):
#   SingleCellExperiment, ggplot2, ggbeeswarm, patchwork, ggh4x, scales, grid
# Optional (only when use_richtext = TRUE): ggtext

plotAbundancesDiffStrip <- function(
    x,
    k                  = "meta20",
    by                 = c("sample_id", "cluster_id"),
    group_by           = "condition",
    group_levels       = NULL,
    n_cols             = 4,
    point_size         = 2,
    facet_ratio        = 1.5,   # aspect ratio (height / width) of each data plot panel
    panel_width_cm     = 3,     # fixed panel WIDTH (cm) across all clusters; locks
                                # data + strip panels to identical width so they align
                                # regardless of how many brackets each cluster has.
    clusters_order     = NULL,
    fill_palette       = NULL,
    # Stats
    external_stats     = NULL,
    show_stats         = TRUE,
    hide_ns            = TRUE,
    # Asterisk / label styling (mirrors add_pvalue_npc convention)
    label_col          = "p.adj.signif",
    label_size         = 3,
    label_fontface     = "plain",
    use_richtext       = TRUE,
    asterisk_pt_multiplier  = 1.4,
    # Asterisk vertical positioning -- scale-invariant. `asterisk_vjust` is in
    # TEXT-HEIGHT units (not data units), so it stays at the same visual offset
    # relative to the bracket no matter how big the figure is.
    #   vjust 0    = text bbox bottom at bracket (text sits just above bar)
    #   vjust 0.5  = text centered on bracket (overlaps bar by half)
    #   vjust 1    = text bbox top at bracket (text below bar)
    # Increase asterisk_vjust to push asterisks DOWN (closer to / overlapping the bar).
    asterisk_vjust     = 0.3,
    label_vjust        = 0,              # vjust for non-asterisk labels (e.g. "ns")
    # DEPRECATED: asterisk_y_offset shifted asterisks in DATA units, which
    # caused fig-height-dependent drift. Kept at 0 default (no effect). If you
    # set it explicitly, it adds on top of the vjust shift and remains
    # scale-dependent. Use asterisk_vjust instead for scale-invariance.
    asterisk_y_offset  = 0,
    bracket_color      = "black",
    bracket_size       = 0.4,
    bracket_tip_npc    = 0.20,           # tip length in strip-local y units (downward)
    label_nudge_npc    = 0.05,           # label vertical offset above the bracket (strip-local y units)
    # Strip layout (NEW: uniform strip across all panels)
    # Every cluster gets a strip of identical height (a fixed fraction of the
    # composite). Within each strip the bracket rows are sized by the cluster
    # with the MOST surviving brackets (`max_brackets`), so spacing between
    # bars is identical across all panels. Clusters with fewer bars get their
    # bars at the BOTTOM of the strip (closer to the data plot below).
    strip_height_fraction = 0.25,        # fraction of composite height for the strip (0.25 = 1/4)
    strip_height_per_bracket = NULL,     # DEPRECATED: variable strip heights -- ignored now;
                                         # kept for back-compat. Use strip_height_fraction instead.
    panel_gap_pt       = 8,              # gap between cluster composites in the grid
    # Y-axis appearance
    y_expand_low_mult  = 0.02,           # bottom y-axis expansion so data isn't pinned to the x-axis
    y_expand_high_mult = 0.05,           # top y-axis expansion (used only when tight_top_axis = FALSE)
    tight_top_axis           = TRUE,     # terminate the y-axis at the "next pretty tick" above the data
                                         # max (with a tiny overhang). Eliminates the dead-space gap
                                         # between the topmost tick and the axis terminus.
    tight_top_axis_overhang  = 0.02,     # multiplicative overhang past the topmost tick (axis_top = top_tick * (1 + this))
    # Title placement
    title_above_strip  = TRUE,           # put the cluster title above the stats strip rather than between strip and data
    # Per-panel typography (patchwork-based plots don't share theme, so they
    # need explicit sizes here rather than letting `+ theme(text=...)` propagate)
    title_size         = 10,             # pt for the cluster title (plot.title)
    axis_text_size     = 8,              # pt for x and y axis tick labels
    # Theme
    nature_style       = TRUE,
    legend_position    = "right",
    # Diagnostic
    verbose            = TRUE          # print estimated min fig.width on call (helps avoid PDF clipping)
) {
    # ---- Package dependencies --------------------------------------------
    if (!requireNamespace("patchwork", quietly = TRUE))
        stop("patchwork is required: install.packages('patchwork').")
    if (!requireNamespace("ggh4x", quietly = TRUE))
        stop("ggh4x is required (force_panelsizes): install.packages('ggh4x').")
    if (isTRUE(use_richtext) && !requireNamespace("ggtext", quietly = TRUE))
        stop("ggtext is required when use_richtext = TRUE: install.packages('ggtext').")

    # ---- Input validation -------------------------------------------------
    by <- match.arg(by)
    stopifnot(
        "panel_width_cm must be a positive number" =
            is.numeric(panel_width_cm) && length(panel_width_cm) == 1 && panel_width_cm > 0,
        "n_cols must be a positive integer" =
            is.numeric(n_cols) && length(n_cols) == 1 && n_cols >= 1,
        "facet_ratio must be a positive number" =
            is.numeric(facet_ratio) && length(facet_ratio) == 1 && facet_ratio > 0,
        "strip_height_fraction must be in (0, 1)" =
            is.numeric(strip_height_fraction) && length(strip_height_fraction) == 1 &&
            strip_height_fraction > 0 && strip_height_fraction < 1
    )

    .wl_check_sce(x, TRUE)

    # ---- Cluster ids ------------------------------------------------------
    # When k matches a colData column, take cluster IDs from there. Otherwise
    # treat k as a cluster_codes resolution name and resolve via the CATALYST-free
    # .wl_check_k() / .wl_cluster_ids() internals (R/wl_internals.R).
    # NB: local var named cluster_id_vec to avoid shadowing .wl_cluster_ids().
    if (k %in% names(colData(x))) {
        cluster_id_vec <- x[[k]]
    } else {
        k <- .wl_check_k(x, k)
        cluster_id_vec <- .wl_cluster_ids(x, k)
    }

    # ---- Build proportion data --------------------------------------------
    ns <- table(cluster_id = cluster_id_vec, sample_id = .wl_sample_ids(x))
    fq <- prop.table(ns, 2) * 100
    df <- as.data.frame(fq)
    m  <- match(df$sample_id, x$sample_id)
    df[[group_by]] <- x[[group_by]][m]

    if (!is.null(clusters_order)) {
        df$cluster_id <- factor(df$cluster_id, levels = clusters_order)
    } else {
        df$cluster_id <- factor(df$cluster_id, levels = unique(as.character(df$cluster_id)))
    }
    clusters <- levels(droplevels(df$cluster_id))

    # ---- Group levels on x-axis -------------------------------------------
    if (is.null(group_levels)) {
        group_levels <- if (is.factor(df[[group_by]])) {
            levels(droplevels(df[[group_by]]))
        } else {
            sort(unique(as.character(df[[group_by]])))
        }
    }
    df[[group_by]] <- factor(df[[group_by]], levels = group_levels)

    # ---- Label transform (asterisk styling) -------------------------------
    # We size up asterisks via a ggtext span. vertical-align CSS is omitted
    # because gridtext only applies it relative to surrounding text on the
    # same line -- and our labels are JUST the asterisk span (no neighbors),
    # so vertical-align had no effect. Vertical positioning of asterisks is
    # done by shifting the data-frame y value directly via `asterisk_y_offset`.
    if (isTRUE(use_richtext)) {
        asterisk_pt <- max(1, round(label_size * 3.5 * asterisk_pt_multiplier))
        label_transform <- function(s) {
            has_ast <- grepl("[*]", s)
            if (any(has_ast)) {
                s_esc <- gsub("\\*", "&#42;", s)
                s[has_ast] <- paste0(
                    "<span style='font-size:", asterisk_pt, "pt; line-height:0.6'>",
                    s_esc[has_ast], "</span>"
                )
            }
            s
        }
    } else {
        label_transform <- function(s) s
    }

    # ---- Builder helpers --------------------------------------------------
    # `show_title` controls whether THIS plot owns the cluster title. When the
    # cluster has a stats strip and `title_above_strip = TRUE`, the strip plot
    # owns the title and the data plot's `show_title` is FALSE.
    build_data_plot <- function(df_cl, cluster_label, show_title) {
        p <- ggplot2::ggplot(
                df_cl,
                ggplot2::aes(x    = .data[[group_by]],
                             y    = Freq,
                             fill = .data[[group_by]])
            ) +
            ggplot2::geom_boxplot(
                color = "grey16", linewidth = 0.5,         # `size` was renamed to `linewidth` in ggplot2 3.4
                alpha = 0.8, outlier.color = NA,
                position = ggplot2::position_dodge(),
                show.legend = TRUE
            ) +
            ggbeeswarm::geom_quasirandom(
                shape = 21, fill = "grey84", size = point_size,
                width = 0.2, alpha = 0.8
            ) +
            ggplot2::labs(x = NULL, y = NULL,
                          title = if (isTRUE(show_title)) cluster_label else NULL) +
            # Y axis: tight-top mode picks the next pretty tick above data max,
            # sets the axis to terminate just past it (no dead space).
            (if (isTRUE(tight_top_axis)) {
                y_max <- max(df_cl$Freq, na.rm = TRUE)
                if (!is.finite(y_max) || y_max <= 0) y_max <- 1
                # pretty() with a slight margin guarantees the top break is
                # strictly above y_max (so "next tick beyond" is well-defined
                # even when y_max already sits on a round number).
                breaks_vec <- pretty(c(0, y_max * 1.001), n = 5)
                top_break  <- max(breaks_vec)
                axis_top   <- top_break * (1 + tight_top_axis_overhang)
                ggplot2::scale_y_continuous(
                    breaks = breaks_vec,
                    limits = c(0, axis_top),
                    expand = ggplot2::expansion(mult = c(y_expand_low_mult, 0))
                )
            } else {
                ggplot2::scale_y_continuous(
                    limits = c(0, NA),
                    expand = ggplot2::expansion(mult = c(y_expand_low_mult, y_expand_high_mult))
                )
            }) +
            ggplot2::scale_x_discrete(limits = group_levels, drop = FALSE) +
            # Lock panel WIDTH only (cm absolute). Height is determined by
            # `theme(aspect.ratio = facet_ratio)` below -- ggh4x's docs say
            # aspect.ratio is honored when only cols is supplied. Setting BOTH
            # cols and rows here caused patchwork to error inside its layout
            # code ("attempt to set an attribute on NULL" / "index out of
            # bounds (unit subsetting)") because rows-in-cm conflict with
            # patchwork's relative-height composition system.
            ggh4x::force_panelsizes(cols = grid::unit(panel_width_cm, "cm")) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                panel.grid       = ggplot2::element_blank(),
                aspect.ratio     = facet_ratio,
                plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5,
                                                         size = title_size,
                                                         margin = ggplot2::margin(b = 2)),
                axis.text        = ggplot2::element_text(color = "black", size = axis_text_size),
                axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1,
                                                         size = axis_text_size),
                axis.ticks       = ggplot2::element_line(color = "black"),
                plot.margin      = ggplot2::margin(2, 2, 2, 2)
                # legend.position is set globally via & below; setting it here
                # too would just be redundant (and panel.spacing is moot in
                # patchwork sub-plots since they have no facets).
            )

        if (!is.null(fill_palette)) {
            p <- p + ggplot2::scale_fill_manual(values = fill_palette, name = group_by, drop = FALSE)
        }

        if (isTRUE(nature_style)) {
            p <- p + ggplot2::theme(
                panel.border     = ggplot2::element_blank(),
                panel.background = ggplot2::element_blank(),
                axis.line        = ggplot2::element_line(color = "black", linewidth = 0.4)
            )
        } else {
            p <- p + ggplot2::theme(
                panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
            )
        }
        p
    }

    build_strip_plot <- function(stat_cl, cluster_label, show_title,
                                 strip_aspect, max_brackets) {
        n_brackets <- if (is.null(stat_cl)) 0 else nrow(stat_cl)
        n_groups   <- length(group_levels)

        # Strip y-axis spans 0..max_brackets across ALL panels. That gives
        # constant bracket spacing within each strip. Bracket rank (1 = top of
        # THIS panel's stack) lives at y = n_brackets - rank + 0.5, so panels
        # with fewer brackets place them at the BOTTOM of the strip (rows
        # nearest the data plot below it).
        y_axis_max <- max(max_brackets, 1)   # avoid degenerate 0 range if no panel has any stats

        # ---- Build the data frames for brackets + labels (may be empty) ----
        if (n_brackets > 0) {
            bracket_y <- n_brackets - seq_len(n_brackets) + 0.5
            x1 <- match(as.character(stat_cl$group1), group_levels)
            x2 <- match(as.character(stat_cl$group2), group_levels)
            x_mid <- (x1 + x2) / 2
            raw_labels <- as.character(stat_cl[[label_col]])
            labels     <- label_transform(raw_labels)

            tip_dy   <- bracket_tip_npc
            nudge_dy <- label_nudge_npc

            seg_h <- data.frame(x = x1, xend = x2, y = bracket_y, yend = bracket_y)
            seg_l <- data.frame(x = x1, xend = x1, y = bracket_y, yend = bracket_y - tip_dy)
            seg_r <- data.frame(x = x2, xend = x2, y = bracket_y, yend = bracket_y - tip_dy)
            seg_df <- rbind(seg_h, seg_l, seg_r)

            is_ast <- grepl("[*]", raw_labels)
            label_y <- bracket_y + nudge_dy
            # `asterisk_y_offset` is a deprecated data-unit shift (default 0
            # = no effect). The scale-invariant control is `asterisk_vjust`
            # below, applied as an aesthetic to the text layer.
            label_y[is_ast] <- label_y[is_ast] + asterisk_y_offset
            vj_vec <- ifelse(is_ast, asterisk_vjust, label_vjust)
            text_df <- data.frame(x = x_mid, y = label_y, label = labels,
                                  vj = vj_vec)
        } else {
            seg_df  <- data.frame(x = numeric(0), xend = numeric(0),
                                  y = numeric(0), yend = numeric(0))
            text_df <- data.frame(x = numeric(0), y = numeric(0),
                                  label = character(0), vj = numeric(0))
        }

        # ---- Plot --------------------------------------------------------
        p <- ggplot2::ggplot()
        if (nrow(seg_df) > 0) {
            p <- p + ggplot2::geom_segment(
                data    = seg_df,
                mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                color   = bracket_color, linewidth = bracket_size
            )
            # vjust is mapped per-row via the `vj` column: asterisks use
            # `asterisk_vjust`, non-asterisks use `label_vjust`. Because vjust
            # is in text-height units (NOT data units), the visual shift is
            # constant in pixels regardless of fig.height -> scale-invariant.
            if (isTRUE(use_richtext)) {
                p <- p + ggtext::geom_richtext(
                    data    = text_df,
                    mapping = ggplot2::aes(x = x, y = y, label = label, vjust = vj),
                    size    = label_size, fontface = label_fontface,
                    label.size    = NA,
                    fill          = NA,
                    label.padding = grid::unit(c(0, 0, 0, 0), "lines")
                )
            } else {
                p <- p + ggplot2::geom_text(
                    data    = text_df,
                    mapping = ggplot2::aes(x = x, y = y, label = label, vjust = vj),
                    size    = label_size, fontface = label_fontface
                )
            }
        }

        p <- p +
            ggplot2::labs(title = if (isTRUE(show_title)) cluster_label else NULL) +
            ggplot2::scale_x_continuous(limits = c(0.5, n_groups + 0.5),
                                        expand = c(0, 0)) +
            # `oob_keep` keeps labels with y < 0 (from large negative
            # asterisk_y_offset) so they still render in the strip-data gutter.
            ggplot2::scale_y_continuous(
                limits = c(0, y_axis_max + 0.5),
                expand = c(0, 0),
                oob    = scales::oob_keep
            ) +
            ggplot2::coord_cartesian(clip = "off") +
            # Lock width only; height comes from `theme(aspect.ratio = strip_aspect)`
            # below. Setting cm rows here too caused patchwork composition to
            # error -- see matching comment in build_data_plot above.
            ggh4x::force_panelsizes(cols = grid::unit(panel_width_cm, "cm")) +
            ggplot2::theme_void() +
            ggplot2::theme(
                aspect.ratio = strip_aspect,
                plot.title   = ggplot2::element_text(face = "bold", hjust = 0.5,
                                                     size = title_size,
                                                     margin = ggplot2::margin(b = 2)),
                plot.margin  = ggplot2::margin(0, 2, 0, 2)
            )
        p
    }

    # ---- Pre-scan: per-cluster filtered stats + global max_brackets -----
    # Determining max_brackets across clusters fixes the strip y-axis range
    # for every panel, which is what gives constant bracket spacing.
    cluster_stats <- lapply(clusters, function(cl) {
        if (!isTRUE(show_stats) || is.null(external_stats)) return(NULL)
        s <- external_stats[as.character(external_stats$cluster_id) == cl, , drop = FALSE]
        if (isTRUE(hide_ns)) {
            s <- s[as.character(s[[label_col]]) != "ns", , drop = FALSE]
        }
        s
    })
    names(cluster_stats) <- clusters
    max_brackets <- max(0L, vapply(cluster_stats,
                                   function(s) if (is.null(s)) 0L else nrow(s),
                                   integer(1)))

    # ---- Strip and data height ratios (uniform across clusters) ----------
    # Every composite is: [strip = strip_height_fraction] + [data = 1 - strip_fraction].
    # Combined with force_panelsizes locking panel widths, this means EVERY
    # cluster's data panel has identical height too.
    strip_h <- strip_height_fraction
    data_h  <- 1 - strip_height_fraction
    strip_aspect <- facet_ratio * strip_h / data_h

    # ---- Build one composite per cluster ----------------------------------
    composite_plots <- lapply(clusters, function(cl) {
        df_cl   <- df[as.character(df$cluster_id) == cl, , drop = FALSE]
        stat_cl <- cluster_stats[[cl]]

        # Strip always exists now; it always owns the cluster title (so titles
        # all sit at the same height regardless of bracket count). The
        # `title_above_strip` arg is retained for back-compat but is effectively
        # the only mode that makes sense with uniform strips.
        data_show_title  <- !isTRUE(title_above_strip)
        strip_show_title <-  isTRUE(title_above_strip)

        p_data  <- build_data_plot(df_cl, cl, show_title = data_show_title)
        p_strip <- build_strip_plot(stat_cl, cluster_label = cl,
                                    show_title   = strip_show_title,
                                    strip_aspect = strip_aspect,
                                    max_brackets = max_brackets)

        # Stack with relative heights via patchwork. Each plot's panel still
        # has a locked width (force_panelsizes cols) and an aspect.ratio
        # constraint on height; plot_layout heights here are the relative
        # weights patchwork uses when allocating cell vertical space. If
        # fig.height is too tight, the panels are constrained by their
        # allocated row height and aspect.ratio may visually saturate.
        p_strip / p_data + patchwork::plot_layout(heights = c(strip_h, data_h))
    })

    # ---- Combine all clusters via wrap_plots ------------------------------
    # `panel_gap_pt` becomes plot.margin on every sub-plot (applied via `&`),
    # creating visible spacing between cluster composites in the grid.
    gap <- panel_gap_pt
    out <- patchwork::wrap_plots(composite_plots, ncol = n_cols) &
        ggplot2::theme(
            legend.position = legend_position,
            plot.margin     = ggplot2::margin(t = gap, r = gap, b = gap, l = gap, unit = "pt")
        )
    if (legend_position != "none") {
        out <- out + patchwork::plot_layout(guides = "collect")
    }

    # ---- Verbose: estimated min fig.width for un-clipped rendering --------
    # force_panelsizes locks panel widths to cm absolute, so if the requested
    # fig.width is too narrow, axes/legends get clipped. Print an estimate so
    # the user can size their f2() / fig.width accordingly. Heuristic only --
    # actual minimum depends on tick-label widths, but this gets you close.
    if (isTRUE(verbose)) {
        axis_width_cm   <- 1.0
        legend_width_cm <- if (legend_position == "right") 3.0 else 0
        gap_cm          <- panel_gap_pt * 0.03528           # pt -> cm
        col_width_cm    <- panel_width_cm + axis_width_cm + 2 * gap_cm
        total_cm        <- n_cols * col_width_cm + legend_width_cm
        total_in        <- total_cm / 2.54
        message(sprintf(
            "plotAbundancesDiffStrip: estimated min fig.width = %.1f in (%.1f cm) [n_cols=%d, panel_width_cm=%g, legend='%s']. Increase if PDFs look clipped.",
            total_in, total_cm, n_cols, panel_width_cm, legend_position
        ))
    }

    out
}
