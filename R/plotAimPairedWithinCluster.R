# plotAimPairedWithinCluster
#
# Self-contained plotter. Plots the proportion of AIM-positive cells
# WITHIN each pheno cluster (e.g. "% of CD4 CTL cells that are AIM-
# responsive") across CMV-condition groups, with unpaired Wilcoxon
# brackets between CMV pairs.
#
# Layout:
#   - facet_wrap(~ cluster_id), one panel per pheno cluster
#   - x = CMV condition, fill = CMV condition
#   - boxplot + jittered donor dots
#   - per-cluster pairwise Wilcoxon brackets above each panel
#   - one figure-level y-axis label on the left (standard ggplot `labs(y = ...)`)
#
# Per-donor math (the "within-cluster" change vs plotAimPaired's default):
#     n_aim   = cells with `act_col == aim_level` in donor × cluster
#     n_tot   = cells in donor × cluster (all act levels)
#     freq    = n_aim / n_tot * 100
#
# Peptide filter:
#   `peptide_filter` is REQUIRED and must be a single peptide level
#   (e.g. "pp65"). Data is filtered to that peptide BEFORE the
#   within-cluster proportions are computed, so the brackets compare CMV
#   groups under one stim condition and the denominator is the pheno
#   cluster count under that same condition.
#
# Returns a ggplot. Pass to f2()/saveFig() to write to disk.

plotAimPairedWithinCluster <- function(
    x,
    k                       = "merging_pheno",
    act_col                 = "merging_act",
    aim_level               = "AIM",
    group_col               = "infectious_status2",
    group_levels            = NULL,
    peptide_col             = "AIM_cond",
    peptide_filter,
    clusters_order          = NULL,
    excluded_clusters       = NULL,
    keep_clusters           = NULL,
    # Stats
    show_stats              = TRUE,
    hide_ns                 = TRUE,
    p_adjust_method         = c("holm", "BH", "bonferroni", "hochberg", "hommel", "BY"),
    # Layout
    n_cols                  = 4,
    point_size              = 1.5,
    log                     = FALSE,
    # Lower-bound floor applied to Freq when log = TRUE. Zero values
    # (donors with no AIM cells in a cluster) → log10(0) = -Inf → ggplot
    # drops them → some (facet × x-group) cells end up with zero plottable
    # points → geom_quasirandom can't compute its range and errors. Floor
    # to `miny` keeps them visible at the bottom of the log axis.
    miny                    = 0.01,
    # Brackets / labels
    label_size              = 3,
    use_richtext            = TRUE,
    asterisk_pt_multiplier  = 1.4,
    asterisk_vjust          = 0.3,
    bracket_size            = 0.4,
    bracket_tip_length      = 0.01,
    bracket_step_increase   = 0.10,    # vertical step between stacked brackets in each panel
    # Typography
    title_size              = 14,      # facet strip title size
    axis_text_size          = 12,
    y_axis_title            = NULL,    # single figure-level y label
    y_axis_title_size       = 14,
    # Palette
    fill_palette            = NULL,
    # Theme
    nature_style            = TRUE,
    legend_position         = "right"
) {
    p_adjust_method <- match.arg(p_adjust_method)

    # ---- Validate ----------------------------------------------------------
    if (missing(peptide_filter) || length(peptide_filter) != 1L) {
        stop("plotAimPairedWithinCluster: `peptide_filter` is required and ",
             "must be a single peptide level (e.g. \"pp65\").")
    }

    cd <- as.data.frame(SingleCellExperiment::colData(x))
    needed <- c("sample_id", k, act_col, group_col, peptide_col)
    missing_cols <- setdiff(needed, names(cd))
    if (length(missing_cols) > 0) {
        stop("Missing required colData columns: ",
             paste(missing_cols, collapse = ", "))
    }
    if (!peptide_filter %in% unique(as.character(cd[[peptide_col]]))) {
        stop("peptide_filter '", peptide_filter, "' not in ", peptide_col, ".")
    }

    # ---- Filter to single peptide -----------------------------------------
    cd <- cd[as.character(cd[[peptide_col]]) == peptide_filter, , drop = FALSE]
    if (nrow(cd) == 0L) stop("After peptide_filter, no cells remain.")

    # ---- Resolve cluster ordering -----------------------------------------
    cluster_factor_levels <- if (is.factor(cd[[k]])) levels(cd[[k]]) else NULL
    all_clusters <- unique(as.character(cd[[k]]))
    if (!is.null(excluded_clusters))
        all_clusters <- setdiff(all_clusters, excluded_clusters)
    if (!is.null(keep_clusters)) {
        missing_keep <- setdiff(keep_clusters, all_clusters)
        if (length(missing_keep) > 0)
            warning("keep_clusters not in data: ",
                    paste(missing_keep, collapse = ", "))
        all_clusters <- intersect(keep_clusters, all_clusters)
    }
    cluster_order_resolved <- if (!is.null(clusters_order)) {
        intersect(clusters_order, all_clusters)
    } else if (!is.null(keep_clusters)) {
        intersect(keep_clusters, all_clusters)
    } else if (!is.null(cluster_factor_levels)) {
        intersect(cluster_factor_levels, all_clusters)
    } else {
        sort(all_clusters)
    }

    # ---- Resolve group levels ---------------------------------------------
    if (is.null(group_levels)) {
        gl_vec <- cd[[group_col]]
        group_levels <- if (is.factor(gl_vec)) levels(droplevels(gl_vec))
                        else sort(unique(as.character(gl_vec)))
    }

    # ---- Per-(sample, cluster) within-cluster AIM% ------------------------
    df <- cd |>
        dplyr::group_by(sample_id, cluster_id = .data[[k]]) |>
        dplyr::summarise(
            n_tot = dplyr::n(),
            n_aim = sum(as.character(.data[[act_col]]) == aim_level),
            Freq  = ifelse(n_tot > 0, n_aim / n_tot * 100, 0),
            .groups = "drop"
        )
    # SEPARATE the raw frequency (used by the Wilcoxon stats) from the
    # plot-side frequency (with optional pseudocount). Stats MUST use the
    # raw value so any future `miny` change — including a per-cluster
    # `auto_miny` mode like plotAimPaired's lines 1138-1156 — can't perturb
    # p-values. The unpaired Wilcoxon is rank-based and a uniform additive
    # shift would preserve ranks today, but a per-cluster shift wouldn't —
    # safer to decouple now than discover the bug later.
    df$Freq_plot <- if (isTRUE(log)) df$Freq + miny else df$Freq

    # Attach group_col (SAFETY: enforce 1:1 sample_id ↔ group_col before
    # the join. If a sample_id appears with multiple group_col values in
    # cd — metadata corruption, e.g. same donor labelled both "Latent" and
    # "Primary" in different rows — `distinct()` yields multiple rows per
    # sample_id and the left_join silently multiplies df rows, duplicating
    # donors in downstream stats. Hard-stop is cheaper than wrong results.
    sample_meta <- cd |>
        dplyr::select(sample_id, dplyr::all_of(group_col)) |>
        dplyr::distinct()
    if (anyDuplicated(sample_meta$sample_id) > 0) {
        dupes <- sample_meta$sample_id[duplicated(sample_meta$sample_id)]
        offenders <- sample_meta[sample_meta$sample_id %in% dupes, , drop = FALSE]
        stop("plotAimPairedWithinCluster: sample_id ↔ ", group_col,
             " is not 1:1. The following sample_ids appear with multiple ",
             group_col, " values in colData:\n",
             paste(utils::capture.output(print(offenders, row.names = FALSE)),
                   collapse = "\n"),
             "\nFix the metadata (one ", group_col, " per sample_id) before plotting.")
    }
    df <- df |> dplyr::left_join(sample_meta, by = "sample_id")
    # Restrict + order. `ordered = FALSE` is explicit because the source
    # column (`merging_pheno` etc.) is often an ORDERED factor, and
    # factor()'s default (`ordered = is.ordered(x)`) would propagate that.
    # Downstream we left_join() this with a stats df that has a PLAIN
    # factor cluster_id, and dplyr refuses to mix ordered + unordered
    # factors. Coercing both sides to plain factors here prevents that.
    df <- df[as.character(df$cluster_id) %in% cluster_order_resolved &
             as.character(df[[group_col]]) %in% group_levels, , drop = FALSE]
    df$cluster_id   <- factor(as.character(df$cluster_id),
                              levels = cluster_order_resolved, ordered = FALSE)
    df[[group_col]] <- factor(as.character(df[[group_col]]),
                              levels = group_levels, ordered = FALSE)

    # Print per-condition sample counts (matches plotAbundancesDiffStrip's
    # convention).
    n_per_cond_df <- df |>
        dplyr::distinct(sample_id, !!rlang::sym(group_col)) |>
        dplyr::count(!!rlang::sym(group_col), name = "n_samples") |>
        as.data.frame()
    message("plotAimPairedWithinCluster: n samples per ", group_col, ":")
    print(n_per_cond_df, row.names = FALSE)

    # ---- Stats: pairwise unpaired Wilcoxon per cluster --------------------
    stats_df <- NULL
    if (isTRUE(show_stats) && length(group_levels) >= 2) {
        pair_mat <- utils::combn(group_levels, 2)
        rows <- list()
        for (cl in cluster_order_resolved) {
            sub_cl <- df[df$cluster_id == cl, , drop = FALSE]
            for (i in seq_len(ncol(pair_mat))) {
                g1 <- pair_mat[1, i]; g2 <- pair_mat[2, i]
                v1 <- sub_cl$Freq[as.character(sub_cl[[group_col]]) == g1]
                v2 <- sub_cl$Freq[as.character(sub_cl[[group_col]]) == g2]
                v1 <- v1[!is.na(v1)]; v2 <- v2[!is.na(v2)]
                if (length(v1) < 2 || length(v2) < 2) {
                    p_val <- NA_real_
                } else {
                    wt <- tryCatch(
                        stats::wilcox.test(v1, v2, paired = FALSE,
                                           exact = FALSE),
                        error = function(e) NULL
                    )
                    p_val <- if (is.null(wt)) NA_real_ else wt$p.value
                }
                rows[[length(rows) + 1L]] <- data.frame(
                    cluster_id = cl,
                    group1     = g1,
                    group2     = g2,
                    n1         = length(v1),
                    n2         = length(v2),
                    p          = p_val,
                    stringsAsFactors = FALSE
                )
            }
        }
        stats_df <- do.call(rbind, rows)
        # Per-cluster correction (rstatix convention)
        stats_df <- stats_df |>
            dplyr::group_by(cluster_id) |>
            dplyr::mutate(
                p.adj = stats::p.adjust(p, method = p_adjust_method)
            ) |>
            dplyr::ungroup()
        # Significance label column
        sig_label <- function(p) {
            ifelse(is.na(p),    "ns",
            ifelse(p < 0.0001, "****",
            ifelse(p < 0.001,  "***",
            ifelse(p < 0.01,   "**",
            ifelse(p < 0.05,   "*",   "ns")))))
        }
        stats_df$p.signif      <- sig_label(stats_df$p)
        stats_df$p.adj.signif  <- sig_label(stats_df$p.adj)
        # Filter ns if requested
        if (isTRUE(hide_ns))
            stats_df <- stats_df[stats_df$p.adj.signif != "ns", , drop = FALSE]
        # Compute y.position per cluster (stacked above each cluster's data
        # max). Needed for stat_pvalue_manual to know where to draw. The
        # positioning math differs for log vs linear y-axis:
        #   linear: y.position = y_max * (1 + step * rank)
        #   log:    y.position = y_max * 10^(step * rank)
        # so visual spacing between stacked brackets stays consistent.
        # Uses Freq_plot (pseudocounted when log = TRUE) so the bracket
        # positions sit above the actual plotted data, not above the raw Freq.
        ymax_per_cluster <- df |>
            dplyr::group_by(cluster_id) |>
            dplyr::summarise(y_max = max(Freq_plot, na.rm = TRUE), .groups = "drop") |>
            dplyr::mutate(
                y_max = ifelse(is.finite(y_max) & y_max > 0,
                               y_max,
                               # Floor so log10(0) doesn't blow up — pick
                               # a value low enough to be visually
                               # irrelevant on the y-axis.
                               if (isTRUE(log)) 0.01 else 1)
            )
        stats_df$cluster_id <- factor(as.character(stats_df$cluster_id),
                                      levels = cluster_order_resolved,
                                      ordered = FALSE)
        if (nrow(stats_df) > 0) {
            stats_df <- stats_df |>
                dplyr::left_join(ymax_per_cluster, by = "cluster_id") |>
                dplyr::group_by(cluster_id) |>
                dplyr::mutate(
                    .rank      = dplyr::row_number(),
                    y.position = if (isTRUE(log)) {
                        y_max * 10^(bracket_step_increase * .rank)
                    } else {
                        y_max * (1 + bracket_step_increase * .rank)
                    }
                ) |>
                dplyr::ungroup() |>
                dplyr::select(-.rank)
        }

        # ---- DIAGNOSTIC: print y.position values so we can SEE if they go
        # ---- absurd (the symptom of the 1e+22 axis blowup). If the numbers
        # ---- printed below look bounded but the rendered axes still extend
        # ---- to 1e+30, the bug is in the SCALE construction, not stats.
        message("plotAimPairedWithinCluster: y.position diagnostic ",
                "(should all be finite, near the data max per cluster):")
        if (!is.null(stats_df) && nrow(stats_df) > 0 &&
            "y.position" %in% names(stats_df)) {
            print(as.data.frame(stats_df[, c("cluster_id", "y_max",
                                             "y.position")]),
                  row.names = FALSE)
        } else {
            message("  (no significant brackets to position)")
        }
        # Console summary
        message("plotAimPairedWithinCluster: stats summary (peptide_filter = '",
                peptide_filter, "'):")
        print(as.data.frame(stats_df[, c("cluster_id", "group1", "group2",
                                         "n1", "n2", "p", "p.adj",
                                         "p.adj.signif")]),
              row.names = FALSE)
    }

    # ---- Plot --------------------------------------------------------------
    # y = Freq_plot (pseudocounted when log = TRUE); the Wilcoxon above used
    # the raw Freq column so stats are unaffected by the plot-side shift.
    p <- ggplot2::ggplot(
            df,
            ggplot2::aes(x    = .data[[group_col]],
                         y    = Freq_plot,
                         fill = .data[[group_col]])
        ) +
        ggplot2::geom_boxplot(
            color = "grey20", linewidth = 0.4,
            alpha = 0.8, outlier.color = NA,
            show.legend = TRUE
        ) +
        # `method = "tukey"` instead of the default `"quasirandom"`: the
        # quasirandom method uses density estimation to compute jitter
        # positions and FAILS WITH `density.default(): need at least 2
        # points to select a bandwidth automatically` whenever any
        # (facet × x-group) cell has only 1 donor — common in this plot
        # because Primary has n=3 and after splitting by cluster some
        # cells end up with n=1. Tukey is a deterministic letter-value
        # placement that doesn't need density estimation and visually
        # matches quasirandom closely at small-to-moderate n. See
        # ledger G17.
        ggbeeswarm::geom_quasirandom(
            shape = 21, color = "black", stroke = 0.4,
            size = point_size, fill = "grey84",
            width = 0.2, alpha = 0.85,
            method = "tukey"
        ) +
        # ggh4x::facet_wrap2 with axes = "all" renders the x-axis (ticks +
        # labels) on EVERY panel — not just the bottom row, which is
        # facet_wrap()'s default. Important here because n_cols × n_rows
        # is uneven for typical merging_pheno cluster counts (12 / 8 → 2
        # rows with 8 + 4 panels), and the trailing row's labels read
        # awkwardly when most panels above lack them. See ledger G16.
        ggh4x::facet_wrap2(~ cluster_id, scales = "free_y", ncol = n_cols,
                           axes = "all") +
        ggplot2::labs(x = NULL, y = y_axis_title, fill = group_col) +
        ggplot2::theme_bw(base_size = axis_text_size) +
        ggplot2::theme(
            strip.background = ggplot2::element_blank(),
            strip.text       = ggplot2::element_text(face = "bold",
                                                     size = title_size),
            axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
            axis.title.y     = ggplot2::element_text(size = y_axis_title_size),
            legend.position  = legend_position,
            panel.grid       = ggplot2::element_blank()
        )

    if (isTRUE(nature_style)) {
        p <- p + ggplot2::theme(
            panel.border = ggplot2::element_blank(),
            axis.line    = ggplot2::element_line(color = "black",
                                                 linewidth = 0.4)
        )
    }

    # ---- Per-facet y-axis scales ------------------------------------------
    # The previous implementation used a single `scale_y_continuous(limits =
    # c(miny, NA))` combined with `facet_wrap2(scales = "free_y")`. That NA
    # upper bound + ggh4x's free_y auto-detection consumed any layer's data
    # (including `stat_pvalue_manual`'s y.position values) when determining
    # per-panel range, producing absurd 1e+22 / 1e+30 upper bounds whenever
    # bracket positioning math touched a panel.
    #
    # Fix (mirrors plotAimPaired's per-cluster axis_top approach, lines
    # ~783-811): compute a FINITE `axis_top` per cluster from the cluster's
    # data max (and bracket positions if any), snap to a nice break value,
    # then apply per-facet scales via `ggh4x::facetted_pos_scales(y = ...)`.
    # That gives every panel its own bounded scale_y_continuous, so brackets
    # outside the range get clipped instead of blowing up the axis.
    #
    # The scales list MUST be in cluster_order_resolved order (the order
    # facets appear), because facetted_pos_scales matches by position.
    # Uses Freq_plot so axis_top accounts for the pseudocount when log=TRUE
    # (matches the y aesthetic mapped above).
    cluster_data_max <- df |>
        dplyr::group_by(cluster_id) |>
        dplyr::summarise(.dmax = max(Freq_plot, na.rm = TRUE), .groups = "drop")

    cluster_bracket_max <- if (!is.null(stats_df) && nrow(stats_df) > 0 &&
                               "y.position" %in% names(stats_df)) {
        stats_df |>
            dplyr::group_by(cluster_id) |>
            dplyr::summarise(.bmax = max(y.position, na.rm = TRUE),
                             .groups = "drop")
    } else {
        data.frame(cluster_id = character(0), .bmax = numeric(0))
    }

    # Per-cluster axis_top. Freq is ALWAYS a percentage (n_aim/n_tot * 100,
    # mathematically bounded to [0, 100]) so the per-cluster axis_top must
    # never exceed 100 — even when bracket positions would otherwise push it
    # past. The bracket positioning math (y_max * 10^(step*rank)) can climb
    # above 100 if y_max is high and there are many brackets; in that case
    # we cap and let stat_pvalue_manual clip the highest brackets.
    PERCENT_CAP <- 100
    cluster_top <- cluster_data_max |>
        dplyr::left_join(cluster_bracket_max, by = "cluster_id") |>
        dplyr::mutate(
            .raw_top = pmax(.dmax,
                            ifelse(is.na(.bmax), -Inf, .bmax),
                            na.rm = TRUE),
            .raw_top = ifelse(is.finite(.raw_top) & .raw_top > 0,
                              .raw_top,
                              if (isTRUE(log)) miny * 10 else 1),
            # Headroom + snap to nice break, then HARD-CAP at PERCENT_CAP.
            axis_top = if (isTRUE(log)) {
                # Next nice log break above raw_top * 1.3, but never above
                # PERCENT_CAP. Snapping uses decade boundaries since this is
                # data on [0.01, 100] — at most three log decades.
                snapped <- 10 ^ ceiling(log10(.raw_top * 1.3))
                pmin(snapped, PERCENT_CAP)
            } else {
                # Linear: 15% headroom, capped at PERCENT_CAP
                pmin(.raw_top * 1.15, PERCENT_CAP)
            }
        )

    # Reorder to match facet order
    cluster_top <- cluster_top[match(cluster_order_resolved,
                                     as.character(cluster_top$cluster_id)), ,
                               drop = FALSE]
    # Guard: any cluster missing from `df` (no rows after filter) — fall back
    cluster_top$axis_top[is.na(cluster_top$axis_top)] <-
        if (isTRUE(log)) miny * 100 else 1

    message("plotAimPairedWithinCluster: per-cluster axis_top:")
    print(as.data.frame(cluster_top[, c("cluster_id", ".dmax", ".bmax",
                                        "axis_top")]),
          row.names = FALSE)

    if (!requireNamespace("ggh4x", quietly = TRUE)) {
        stop("ggh4x is required for facetted_pos_scales: ",
             "install.packages('ggh4x').")
    }

    # Plain-number formatter for log-axis breaks: 0.01 / 0.1 / 1 / 10 / 100
    # (ggplot's default formatter renders these as 1e-02 / 1e-01 / 1e+00 /
    # 1e+01 / 1e+02 which is ugly for percentage data). Ported from
    # plotAimPaired.R line 271.
    .fmt_log_break <- function(x) {
        vapply(x, function(v) {
            s <- format(v, scientific = FALSE, trim = TRUE)
            s <- sub("(\\..*?)0+$", "\\1", s)  # ".10" -> ".1"; ".00" -> "."
            s <- sub("\\.$",        "",   s)   # trailing "." -> ""
            s
        }, character(1))
    }

    if (isTRUE(log)) {
        y_scales <- lapply(cluster_top$axis_top, function(top) {
            # Explicit decade breaks within [miny, top] so labels are
            # predictable: e.g. 0.01, 0.1, 1, 10, 100.
            decade_lo  <- floor(log10(miny))
            decade_hi  <- ceiling(log10(top))
            breaks_vec <- 10 ^ seq(decade_lo, decade_hi)
            breaks_vec <- breaks_vec[breaks_vec >= miny & breaks_vec <= top]
            ggplot2::scale_y_continuous(
                trans  = "log10",
                limits = c(miny, top),
                breaks = breaks_vec,
                labels = .fmt_log_break(breaks_vec),
                expand = ggplot2::expansion(mult = c(0, 0))
            )
        })
    } else {
        y_scales <- lapply(cluster_top$axis_top, function(top) {
            ggplot2::scale_y_continuous(
                limits = c(0, top),
                expand = ggplot2::expansion(mult = c(0, 0))
            )
        })
    }

    p <- p + ggh4x::facetted_pos_scales(y = y_scales)

    if (!is.null(fill_palette)) {
        p <- p + ggplot2::scale_fill_manual(values = fill_palette)
    }

    # ---- Brackets via ggpubr::stat_pvalue_manual --------------------------
    # stat_pvalue_manual renders labels via geom_text (NOT geom_richtext),
    # so wrapping the asterisks in `<span style='font-size:Xpt'>...</span>`
    # would print the literal HTML tags. We use plain "*"/"**"/etc. labels
    # and control their pt size directly via `size` (multiplied by
    # asterisk_pt_multiplier when the label is an asterisk row).
    #
    # `inherit.aes = FALSE` is critical: the base ggplot has a global
    # `aes(x = .data[[group_col]], ...)` and stat_pvalue_manual would
    # otherwise try to evaluate that on stats_df (which doesn't have the
    # group_col column), producing "Column `infectious_status2` not found".
    if (!is.null(stats_df) && nrow(stats_df) > 0) {
        if (!requireNamespace("ggpubr", quietly = TRUE))
            stop("ggpubr is required for stats brackets: ",
                 "install.packages('ggpubr').")

        # Apply pt-multiplier to the rendered label size, applied uniformly.
        # Asterisks render at `label_size * asterisk_pt_multiplier`,
        # non-asterisks at `label_size`.
        # NB: ggpubr's `size` is geom_text units (mm), same as label_size's
        # native scale. Multipliers work as expected.
        ast_size <- label_size * asterisk_pt_multiplier
        is_ast_row <- grepl("[*]", stats_df$p.adj.signif)

        if (any(is_ast_row)) {
            p <- p + ggpubr::stat_pvalue_manual(
                stats_df[is_ast_row, , drop = FALSE],
                label        = "p.adj.signif",
                tip.length   = bracket_tip_length,
                bracket.size = bracket_size,
                size         = ast_size,
                vjust        = asterisk_vjust,
                inherit.aes  = FALSE
            )
        }
        if (any(!is_ast_row)) {
            p <- p + ggpubr::stat_pvalue_manual(
                stats_df[!is_ast_row, , drop = FALSE],
                label        = "p.adj.signif",
                tip.length   = bracket_tip_length,
                bracket.size = bracket_size,
                size         = label_size,
                vjust        = asterisk_vjust,
                inherit.aes  = FALSE
            )
        }
    }

    # Stash source data + stats for f2()'s saveExcel.
    attr(p, "source_data")  <- df
    attr(p, "source_stats") <- stats_df
    attr(p, "n_per_condition") <- n_per_cond_df

    p
}
