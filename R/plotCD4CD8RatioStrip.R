# plotCD4CD8RatioStrip
#
# Sibling of plotCD4CD8Ratio with the separated stats-strip layout. Single-
# panel case: a [stats strip / data plot] composite drawn via patchwork.
#
# Mirrors the upstream function in computing per-sample ratio = numerator /
# denominator (or (numerator + pc) / (denominator + pc) if pseudocount > 0),
# but adds:
#   - stats strip above the data plot, with uniform spacing if multiple
#     brackets fit
#   - rounded bracket caps (no exposed corners)
#   - scale-invariant asterisk positioning via asterisk_vjust
#   - tight_top_axis so the y-axis terminates at the next pretty tick above
#     the highest data point (no dead space)
#   - attr(p, "source_data") + attr(p, "source_stats") for f2(saveExcel=TRUE)
#     long-format export (full unfiltered p-values included)
#
# Reference line at y = 1 (unity ratio) kept from the upstream function.
#
# Dependencies (loaded by SETUP 1 of the calling .Rmd):
#   SingleCellExperiment, ggplot2, ggbeeswarm, patchwork, ggh4x, scales, grid,
#   dplyr, tidyr, rstatix. Optional: ggtext (when use_richtext = TRUE).

plotCD4CD8RatioStrip <- function(
    x,
    k = "merging1",
    group_by            = "condition3",
    patient_by          = "sample_id",
    shape_by            = NULL,
    numerator_cluster   = "CD4",
    denominator_cluster = "CD8",
    pseudocount         = 0,
    group_levels        = NULL,
    fill_palette        = NULL,
    point_size          = 2,
    facet_ratio         = 1.5,
    panel_width_cm      = 4,
    # Log axis (mirrors upstream)
    log                 = FALSE,
    miny                = 0.01,
    maxy                = NA,
    # Stats
    external_stats      = NULL,
    show_stats          = TRUE,
    hide_ns             = FALSE,
    # Multiple-testing correction view. Three options; meaningful only when
    # `external_stats` carries the corresponding signif columns. For internally
    # computed stats (external_stats = NULL) there's only one cluster
    # being tested, so per-pair / per-cluster / global all yield the same
    # family — switch is a no-op there.
    #   "per_pair"    -> diffcyt default. Uses p.adj.signif.
    #   "per_cluster" -> recommended (per PI); BH only across the pairs for
    #                   this cluster. Uses p.adj_per_cluster.signif.
    #   "global"      -> BH across the full pre-supplied family.
    #                   Uses p.adj_global.signif.
    correction          = c("per_pair", "per_cluster", "global"),
    # Asterisk / label styling
    label_col           = "p.adj.signif",
    label_size          = 3,
    label_fontface      = "plain",
    use_richtext        = TRUE,
    asterisk_pt_multiplier = 1.4,
    asterisk_vjust      = 0.3,
    label_vjust         = 0,
    asterisk_y_offset   = 0,
    bracket_color       = "black",
    bracket_size        = 0.4,
    bracket_tip_npc     = 0.20,
    label_nudge_npc     = 0.05,
    # Strip layout
    strip_height_fraction = 0.25,
    panel_gap_pt          = 8,
    # Y-axis
    y_expand_low_mult   = 0.02,
    y_expand_high_mult  = 0.05,
    tight_top_axis           = TRUE,
    tight_top_axis_overhang  = 0.02,
    # Title / typography
    title               = NULL,
    title_above_strip   = TRUE,
    title_size          = 12,
    axis_text_size      = 10,
    # Theme
    nature_style        = TRUE,
    legend_position     = "right",
    # Internal stats method (when external_stats is NULL)
    stats_method        = c("dunn", "wilcox"),
    p_adjust_method     = "holm",
    verbose             = TRUE
) {
    # ---- Package dependencies --------------------------------------------
    if (!requireNamespace("patchwork", quietly = TRUE))
        stop("patchwork is required: install.packages('patchwork').")
    if (!requireNamespace("ggh4x", quietly = TRUE))
        stop("ggh4x is required (force_panelsizes): install.packages('ggh4x').")
    if (!requireNamespace("rstatix", quietly = TRUE) && isTRUE(show_stats) && is.null(external_stats))
        stop("rstatix is required when computing stats internally: install.packages('rstatix').")
    if (isTRUE(use_richtext) && !requireNamespace("ggtext", quietly = TRUE))
        stop("ggtext is required when use_richtext = TRUE: install.packages('ggtext').")
    stats_method <- match.arg(stats_method)

    # Correction-view convenience switch (see param doc above). Only meaningful
    # when external_stats carries the global column. For internal stats this
    # is effectively a no-op.
    correction <- match.arg(correction)
    if (identical(label_col, "p.adj.signif")) {
        label_col <- switch(correction,
                            per_pair    = "p.adj.signif",
                            per_cluster = "p.adj_per_cluster.signif",
                            global      = "p.adj_global.signif")
    }

    if (!inherits(x, "SingleCellExperiment"))
        stop("x must be a SingleCellExperiment object.")

    df <- as.data.frame(SingleCellExperiment::colData(x))

    # ---- Resolve k column (handle cluster_codes case) --------------------
    if (!(k %in% colnames(df))) {
        # CATALYST-free resolver (.wl_cluster_ids reads metadata(x)$cluster_codes).
        resolved_k <- tryCatch(
            as.character(.wl_cluster_ids(x, k = k)),
            error = function(e) NULL
        )
        if (is.null(resolved_k)) {
            cluster_codes <- SingleCellExperiment::metadata(x)$cluster_codes
            if (!is.null(cluster_codes) && k %in% colnames(cluster_codes) &&
                "cluster_id" %in% colnames(df)) {
                id_vec <- as.character(df$cluster_id)
                idx    <- match(id_vec, as.character(cluster_codes$cluster_id))
                resolved_k <- as.character(cluster_codes[[k]][idx])
            }
        }
        if (is.null(resolved_k) || all(is.na(resolved_k))) {
            stop("Missing required colData column: ", k,
                 ". Could not derive from .wl_cluster_ids or metadata(x)$cluster_codes.")
        }
        df[[k]] <- resolved_k
    }

    # ---- Validate inputs -------------------------------------------------
    needed <- c(k, group_by, patient_by)
    missing_cols <- setdiff(needed, colnames(df))
    if (length(missing_cols))
        stop("Missing required colData columns: ", paste(missing_cols, collapse = ", "))

    if (!is.null(shape_by) && !(shape_by %in% colnames(df)))
        stop("shape_by column not found in colData: ", shape_by)

    if (identical(numerator_cluster, denominator_cluster))
        stop("numerator_cluster and denominator_cluster must be different.")

    # ---- Filter to the two clusters of interest --------------------------
    df$cluster_lab <- as.character(df[[k]])
    obs <- unique(df$cluster_lab)
    if (!(numerator_cluster   %in% obs)) warning("numerator_cluster not found in ", k, ": ", numerator_cluster)
    if (!(denominator_cluster %in% obs)) warning("denominator_cluster not found in ", k, ": ", denominator_cluster)

    df <- df[df$cluster_lab %in% c(numerator_cluster, denominator_cluster), , drop = FALSE]
    if (nrow(df) == 0)
        stop("No cells found with ", numerator_cluster, " or ", denominator_cluster, " in ", k, ".")

    # ---- Compute ratio per (patient, group) -----------------------------
    group_cols <- c(patient_by, group_by)
    if (!is.null(shape_by)) group_cols <- c(group_cols, shape_by)

    count_df <- df %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "cluster_lab")))) %>%
        dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")

    ratio_df <- count_df %>%
        tidyr::pivot_wider(names_from = cluster_lab, values_from = n_cells, values_fill = 0)

    if (!(numerator_cluster %in% colnames(ratio_df))) ratio_df[[numerator_cluster]] <- 0
    if (!(denominator_cluster %in% colnames(ratio_df))) ratio_df[[denominator_cluster]] <- 0

    num_n <- ratio_df[[numerator_cluster]]
    den_n <- ratio_df[[denominator_cluster]]
    if (pseudocount == 0) {
        ratio_df$ratio <- ifelse(den_n == 0, NA_real_, num_n / den_n)
    } else {
        ratio_df$ratio <- (num_n + pseudocount) / (den_n + pseudocount)
    }
    names(ratio_df)[names(ratio_df) == numerator_cluster]   <- "numerator_n"
    names(ratio_df)[names(ratio_df) == denominator_cluster] <- "denominator_n"

    if (log) {
        # Same convention as upstream: shift by miny so zeros plot at miny.
        ratio_df$ratio <- ratio_df$ratio + miny
    }

    # ---- Group levels for x-axis ----------------------------------------
    if (is.null(group_levels)) {
        group_levels <- if (is.factor(ratio_df[[group_by]])) {
            levels(droplevels(ratio_df[[group_by]]))
        } else {
            sort(unique(as.character(ratio_df[[group_by]])))
        }
    }
    ratio_df[[group_by]] <- factor(ratio_df[[group_by]], levels = group_levels)

    # ---- n samples per condition (Immunity figure legend reporting) -------
    # `patient_by` is the per-sample id used to build ratio_df, so counting
    # distinct patient_by values per group_by gives the per-condition n.
    n_per_cond_df <- ratio_df %>%
        dplyr::distinct(!!rlang::sym(patient_by), !!rlang::sym(group_by)) %>%
        dplyr::count(!!rlang::sym(group_by), name = "n_samples") %>%
        as.data.frame()
    message("plotCD4CD8RatioStrip: n samples per ", group_by, ":")
    print(n_per_cond_df, row.names = FALSE)

    # ---- Stats ----------------------------------------------------------
    full_stats <- NULL
    rendered_stats <- NULL
    if (isTRUE(show_stats)) {
        stat_df <- ratio_df[!is.na(ratio_df$ratio), , drop = FALSE]
        if (!is.null(external_stats)) {
            full_stats <- as.data.frame(external_stats)
        } else if (nrow(stat_df) >= 2) {
            full_stats <- tryCatch({
                if (stats_method == "dunn") {
                    rstatix::dunn_test(stat_df,
                                       as.formula(paste("ratio ~", group_by)),
                                       p.adjust.method = p_adjust_method)
                } else {
                    rstatix::wilcox_test(stat_df,
                                         as.formula(paste("ratio ~", group_by)),
                                         p.adjust.method = p_adjust_method)
                }
            }, error = function(e) NULL)
        }
        if (!is.null(full_stats) && nrow(full_stats) > 0) {
            rendered_stats <- full_stats
            if (isTRUE(hide_ns)) {
                rendered_stats <- rendered_stats[as.character(rendered_stats[[label_col]]) != "ns", , drop = FALSE]
            }
        }
    }
    n_brackets <- if (is.null(rendered_stats)) 0 else nrow(rendered_stats)

    # ---- Strip / data height fractions (single-panel) -------------------
    strip_h <- strip_height_fraction
    data_h  <- 1 - strip_height_fraction
    strip_aspect <- facet_ratio * strip_h / data_h

    # ---- Label transform (asterisk styling) -----------------------------
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

    # ---- Data plot ------------------------------------------------------
    y_label <- paste0(numerator_cluster, "/", denominator_cluster, " ratio",
                      if (log) " (log10)" else "")

    p_data <- ggplot2::ggplot(ratio_df,
                              ggplot2::aes(x    = .data[[group_by]],
                                           y    = ratio,
                                           fill = .data[[group_by]])) +
        ggplot2::geom_hline(yintercept = 1, linetype = "dotted",
                            color = "grey40", linewidth = 0.5) +
        ggplot2::geom_boxplot(
            color = "grey16", linewidth = 0.5, width = 0.75,
            alpha = 0.8, outlier.color = NA, show.legend = TRUE
        )

    if (!is.null(shape_by)) {
        p_data <- p_data + ggbeeswarm::geom_quasirandom(
            mapping = ggplot2::aes(shape = .data[[shape_by]]),
            size = point_size, width = 0.2
        )
    } else {
        p_data <- p_data + ggbeeswarm::geom_quasirandom(
            shape = 21, fill = "grey84", size = point_size,
            width = 0.2, alpha = 0.8
        )
    }

    p_data <- p_data +
        ggplot2::labs(x = NULL, y = y_label,
                      title = if (!isTRUE(title_above_strip)) title else NULL) +
        ggplot2::scale_x_discrete(limits = group_levels, drop = FALSE) +
        ggh4x::force_panelsizes(cols = grid::unit(panel_width_cm, "cm")) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            panel.grid       = ggplot2::element_blank(),
            aspect.ratio     = facet_ratio,
            plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5,
                                                     size = title_size,
                                                     margin = ggplot2::margin(b = 2)),
            axis.title.y     = ggplot2::element_text(size = title_size),
            axis.text        = ggplot2::element_text(color = "black", size = axis_text_size),
            axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1,
                                                     size = axis_text_size),
            axis.ticks       = ggplot2::element_line(color = "black"),
            plot.margin      = ggplot2::margin(2, 2, 2, 2)
        )

    if (log) {
        # Tight-top on a log10 axis: terminate at the next decade above the
        # data max (the natural "next tick" for log axes). Apply the same
        # `tight_top_axis_overhang` multiplier as the linear path for a
        # tiny visible margin past the topmost tick. If `maxy` was set
        # explicitly by the caller, respect that and skip the tight-top.
        decade_breaks <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)
        if (isTRUE(tight_top_axis) && (is.na(maxy) || is.null(maxy))) {
            y_max <- max(ratio_df$ratio, na.rm = TRUE)
            if (!is.finite(y_max) || y_max <= 0) y_max <- 1
            top_decade <- 10 ^ ceiling(log10(y_max * 1.001))
            axis_top   <- top_decade * (1 + tight_top_axis_overhang)
            breaks_vec <- decade_breaks[decade_breaks >= miny &
                                        decade_breaks <= top_decade]
            p_data <- p_data + ggplot2::scale_y_continuous(
                trans  = "log10",
                limits = c(miny, axis_top),
                breaks = breaks_vec,
                labels = breaks_vec
            )
        } else {
            # Caller pinned maxy, or tight_top off -> legacy behavior.
            breaks_vec <- decade_breaks[decade_breaks >= miny &
                                        (is.na(maxy) | decade_breaks <= maxy)]
            p_data <- p_data + ggplot2::scale_y_continuous(
                trans  = "log10",
                limits = c(miny, maxy),
                breaks = breaks_vec,
                labels = breaks_vec
            )
        }
    } else if (isTRUE(tight_top_axis)) {
        y_max <- max(ratio_df$ratio, na.rm = TRUE)
        if (!is.finite(y_max) || y_max <= 0) y_max <- 1
        breaks_vec <- pretty(c(0, y_max * 1.001), n = 5)
        top_break  <- max(breaks_vec)
        axis_top   <- top_break * (1 + tight_top_axis_overhang)
        p_data <- p_data + ggplot2::scale_y_continuous(
            breaks = breaks_vec,
            limits = c(0, axis_top),
            expand = ggplot2::expansion(mult = c(y_expand_low_mult, 0))
        )
    } else {
        p_data <- p_data + ggplot2::scale_y_continuous(
            limits = c(0, NA),
            expand = ggplot2::expansion(mult = c(y_expand_low_mult, y_expand_high_mult))
        )
    }

    if (!is.null(fill_palette)) {
        p_data <- p_data + ggplot2::scale_fill_manual(values = fill_palette, name = group_by, drop = FALSE)
    }

    if (isTRUE(nature_style)) {
        p_data <- p_data + ggplot2::theme(
            panel.border     = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            axis.line        = ggplot2::element_line(color = "black", linewidth = 0.4)
        )
    } else {
        p_data <- p_data + ggplot2::theme(
            panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
        )
    }

    # ---- Strip plot (uniform with the other strip functions) ------------
    n_groups   <- length(group_levels)
    y_axis_max <- max(n_brackets, 1)

    if (n_brackets > 0) {
        bracket_y <- n_brackets - seq_len(n_brackets) + 0.5
        x1 <- match(as.character(rendered_stats$group1), group_levels)
        x2 <- match(as.character(rendered_stats$group2), group_levels)
        x_mid <- (x1 + x2) / 2
        raw_labels <- as.character(rendered_stats[[label_col]])
        labels     <- label_transform(raw_labels)

        seg_h <- data.frame(x = x1, xend = x2, y = bracket_y, yend = bracket_y)
        seg_l <- data.frame(x = x1, xend = x1, y = bracket_y, yend = bracket_y - bracket_tip_npc)
        seg_r <- data.frame(x = x2, xend = x2, y = bracket_y, yend = bracket_y - bracket_tip_npc)
        seg_df <- rbind(seg_h, seg_l, seg_r)

        is_ast <- grepl("[*]", raw_labels)
        label_y <- bracket_y + label_nudge_npc
        label_y[is_ast] <- label_y[is_ast] + asterisk_y_offset
        vj_vec  <- ifelse(is_ast, asterisk_vjust, label_vjust)
        text_df <- data.frame(x = x_mid, y = label_y, label = labels, vj = vj_vec)
    } else {
        seg_df  <- data.frame(x = numeric(0), xend = numeric(0), y = numeric(0), yend = numeric(0))
        text_df <- data.frame(x = numeric(0), y = numeric(0), label = character(0), vj = numeric(0))
    }

    p_strip <- ggplot2::ggplot()
    if (nrow(seg_df) > 0) {
        p_strip <- p_strip + ggplot2::geom_segment(
            data    = seg_df,
            mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
            color   = bracket_color, linewidth = bracket_size,
            lineend = "round", linejoin = "round"
        )
        if (isTRUE(use_richtext)) {
            p_strip <- p_strip + ggtext::geom_richtext(
                data    = text_df,
                mapping = ggplot2::aes(x = x, y = y, label = label, vjust = vj),
                size    = label_size, fontface = label_fontface,
                label.size = NA, fill = NA,
                label.padding = grid::unit(c(0, 0, 0, 0), "lines")
            )
        } else {
            p_strip <- p_strip + ggplot2::geom_text(
                data    = text_df,
                mapping = ggplot2::aes(x = x, y = y, label = label, vjust = vj),
                size    = label_size, fontface = label_fontface
            )
        }
    }

    p_strip <- p_strip +
        ggplot2::labs(title = if (isTRUE(title_above_strip)) title else NULL) +
        ggplot2::scale_x_continuous(limits = c(0.5, n_groups + 0.5), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(limits = c(0, y_axis_max + 0.5),
                                    expand = c(0, 0),
                                    oob    = scales::oob_keep) +
        ggplot2::coord_cartesian(clip = "off") +
        ggh4x::force_panelsizes(cols = grid::unit(panel_width_cm, "cm")) +
        ggplot2::theme_void() +
        ggplot2::theme(
            aspect.ratio = strip_aspect,
            plot.title   = ggplot2::element_text(face = "bold", hjust = 0.5,
                                                 size = title_size,
                                                 margin = ggplot2::margin(b = 2)),
            plot.margin  = ggplot2::margin(0, 2, 0, 2)
        )

    # ---- Compose composite ----------------------------------------------
    gap <- panel_gap_pt
    out <- p_strip / p_data +
        patchwork::plot_layout(heights = c(strip_h, data_h))
    out <- out & ggplot2::theme(
        legend.position = legend_position,
        plot.margin     = ggplot2::margin(t = gap, r = gap, b = gap, l = gap, unit = "pt")
    )
    if (legend_position != "none") {
        out <- out + patchwork::plot_layout(guides = "collect")
    }

    # ---- Stash data + full unfiltered stats for f2() --------------------
    attr(out, "source_data")  <- ratio_df
    attr(out, "n_per_condition") <- n_per_cond_df
    attr(out, "source_stats") <- full_stats

    if (isTRUE(verbose)) {
        axis_width_cm   <- 1.0
        legend_width_cm <- if (legend_position == "right") 3.0 else 0
        gap_cm          <- panel_gap_pt * 0.03528
        total_cm        <- panel_width_cm + axis_width_cm + legend_width_cm + 2 * gap_cm
        total_in        <- total_cm / 2.54
        message(sprintf(
            "plotCD4CD8RatioStrip: estimated min fig.width = %.1f in (%.1f cm) [panel_width_cm=%g, legend='%s']. Increase if PDFs look clipped.",
            total_in, total_cm, panel_width_cm, legend_position
        ))
    }

    out
}
