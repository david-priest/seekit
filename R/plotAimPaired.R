# plotAimPaired
#
# AIM-style paired slope-graph plotter with Wilcoxon signed-rank paired stats.
#
# The paper-standard layout for AIM data (e.g. Grifoni et al. 2020, Cell Fig 2C):
# show raw % AIM-responding cells per donor, with **lines connecting each donor**
# across stimulation conditions (No Peptide -> peptide), and a paired Wilcoxon
# p-value/star above each "reference -> peptide" comparison.
#
# For each cluster (k level) this function produces a row of mini-panels --
# one per group_col level (e.g. CMV status: Never / Latent / Persistent IgM /
# Primary). Within each mini-panel, donors are plotted at the peptide-condition
# x positions (No Peptide / IE1 / pp65) with thin lines connecting the same
# donor's points. Two paired Wilcoxon tests per mini-panel are computed:
# reference vs each non-reference peptide (e.g. No Peptide vs IE1; No Peptide
# vs pp65). The stars sit in a stats strip above the mini-panel, matching the
# existing strip plotters' convention.
#
# All clusters are wrapped via patchwork into the final figure.
#
# Multiple-testing correction (matches the rest of the paper):
#   - p.adj          -> per-cluster BH across the cluster's
#                       (group_col levels x non-reference peptides) tests
#   - p.adj_global   -> BH across every test in the figure
# The `correction` param picks which column the rendered stars come from
# (default "per_cluster"); the unfiltered stats are still saved on the
# returned patchwork as attr(., "source_stats").
#
# Dependencies: SingleCellExperiment, ggplot2, ggbeeswarm, patchwork, ggh4x,
# scales, grid, dplyr, tidyr, rstatix. Optional: ggtext (use_richtext = TRUE).

plotAimPaired <- function(
    x,
    k,                                # colData column with cluster labels (or CATALYST k name)
    group_col       = "infectious_status2",   # mini-panel facet (one per level)
    peptide_col     = "AIM_cond",             # x-axis of each mini-panel
    patient_col     = "patient_id",           # paired-by
    reference_level = "No Peptide",           # the peptide_col level used as the within-donor baseline
    # ---- "WITHIN-CLUSTER" denominator mode -----------------------------------
    # Default behaviour: y = (cells in cluster `k` level) / (total cells in
    # sample) — i.e. the fraction of the parent lineage in that cluster.
    #
    # When `act_col` + `aim_level` are BOTH non-NULL the denominator AND
    # numerator change:
    #   - panels still defined by `k` (e.g. one panel per pheno cluster)
    #   - numerator = cells with `act_col == aim_level` within that panel
    #     (i.e. AIM-positive cells inside the pheno)
    #   - denominator = all cells in that panel (per sample)
    # So y = AIM%-within-cluster, not AIM%-of-lineage.
    #
    # Typical use: k = "merging_pheno", act_col = "merging_act",
    # aim_level = "AIM" → "AIM-positive % per pheno cluster". Consumed by
    # the wrapper plotAimPairedWithinCluster().
    act_col         = NULL,
    aim_level       = NULL,
    group_levels    = NULL,                   # explicit group_col order; defaults to factor levels
    peptide_levels  = NULL,                   # explicit peptide_col order; defaults to factor levels
    clusters_order  = NULL,
    excluded_clusters = NULL,
    # Whitelist filter applied AFTER per-donor proportions are calculated, so
    # denominators (total cells per sample) remain unaffected -- only the SET
    # of clusters whose panels get rendered is reduced. Use case: paper-figure
    # subsets of merging_combined plots (e.g. "show only the three AIM
    # cluster panels that go into Fig 5"). NULL keeps all clusters (current
    # behaviour). If supplied AND `clusters_order` is NULL, the order of
    # `keep_clusters` is used as the cluster panel order; pass `clusters_order`
    # explicitly to override.
    keep_clusters     = NULL,
    merging_col     = NULL,                   # forwarded to plotAbundancesA1-style cluster_id resolution
    # Stats
    show_stats      = TRUE,
    hide_ns         = FALSE,
    correction      = c("per_cluster", "global"),
    label_col       = "p.adj.signif",         # override of `correction` if not the default value
    # Multiple-testing correction method (passed straight to stats::p.adjust).
    # Default "holm" is appropriate for SMALL pre-specified families like
    # ours (per-cluster: 4 CMV statuses x 2 peptide contrasts = 8 tests),
    # provides FWER control rather than BH's FDR control, and matches
    # rstatix::wilcox_test's default for paired-comparison families. Switch
    # to "BH" for consistency with the diffcyt-based ex vivo analyses.
    p_adjust_method = c("holm", "BH", "bonferroni", "hochberg", "hommel", "BY"),
    # Geometry
    n_cols                 = 4,
    panel_width_cm         = 3,
    facet_ratio            = 1.5,
    strip_height_fraction  = 0.25,
    panel_gap_pt           = 8,
    point_size             = 2,
    line_alpha             = 0.4,
    line_color             = "grey50",
    line_size              = 0.3,
    # Boxplots behind the points + lines. Default ON to give a visual sense
    # of distribution; set FALSE for clean Grifoni-style "points + donor lines"
    # only.
    show_boxplots          = TRUE,
    boxplot_alpha          = 0.6,
    boxplot_linewidth      = 0.3,
    boxplot_width          = 0.6,
    # Deterministic within-x-slot jitter so multiple donors don't stack and
    # so the donor lines terminate exactly on the points (the line frame
    # inherits .x_pos from the same precomputed jitter).
    jitter_width           = 0.08,
    jitter_seed            = 42,
    # Y-axis
    y_expand_low_mult        = 0.02,
    y_expand_high_mult       = 0.05,
    tight_top_axis           = TRUE,
    tight_top_axis_overhang  = 0.02,
    log                      = FALSE,
    miny                     = 0.01,
    # When TRUE (default), the y-axis lower limit AND the zero-shift are
    # raised to 0.1 in any cluster whose minimum non-zero frequency is
    # >= 0.1 -- avoids wasted bottom log-space when the data never goes
    # near 0.01. Falls back to the `miny` parameter above when the cluster
    # has data below 0.1. Per-cluster: in shared_y_per_row mode each
    # cluster row can pick its own appropriate bottom.
    auto_miny                = TRUE,
    # On log axes, how should the y-axis terminate above the data max?
    #   "nice_2_5" -> smallest "nice" log break (..., 0.5, 1, 2, 5, 10, 20, 50, ...)
    #                above y_max * (1 + overhang). Avoids huge whitespace when
    #                the data sits just inside the current decade (e.g. y_max
    #                = 1.2 terminates at 2 instead of 10).
    #   "decade"   -> snap to the next power of 10 above y_max (old behaviour).
    log_axis_top_style       = c("nice_2_5", "decade"),
    # Asterisk / bracket styling (mirrors plotAbundancesDiffStrip)
    label_size               = 3,
    label_fontface           = "plain",
    use_richtext             = TRUE,
    asterisk_pt_multiplier   = 1.4,
    asterisk_vjust           = 0.3,
    label_vjust              = 0,
    asterisk_y_offset        = 0,
    bracket_color            = "black",
    bracket_size             = 0.4,
    bracket_tip_npc          = 0.20,
    label_nudge_npc          = 0.05,
    # Title / typography
    title_above_strip = TRUE,
    title_size        = 12,
    axis_text_size    = 10,
    row_label_size    = NULL,                  # defaults to title_size
    # Theme
    fill_palette      = NULL,
    nature_style      = TRUE,
    legend_position   = "right",
    # Legend appearance (right-side legend by default).
    # legend_key_size: legend-glyph diameter in mm (was effectively
    # ~3-4mm from ggplot defaults). legend_text_size /
    # legend_title_size override axis_text_size / row_label_size for the
    # legend specifically; legend_spacing_y_pt tightens vertical
    # whitespace between entries.
    legend_key_size     = 6,
    legend_text_size    = NULL,  # NULL -> falls back to axis_text_size
    legend_title_size   = NULL,  # NULL -> falls back to row_label_size
    legend_spacing_y_pt = 2,
    # Shared-y / tight-row layout. When TRUE all panels within the same
    # cluster row share the same y-axis range (max across that cluster's
    # groups), and y-axis chrome is hidden on every panel except the
    # leftmost one. This matches the look of plotAbundancesA1's facet_wrap
    # output. tighten_within_row reduces the left/right plot.margin so
    # adjacent same-row panels sit close together (vertical row gap is
    # preserved).
    shared_y_per_row  = TRUE,
    tighten_within_row = TRUE,
    panel_gap_x_pt    = 2,
    # Manually-drawn x-axis line that extends past each panel's boundaries,
    # so adjacent panels' lines overlap and read as one continuous axis
    # (facet_wrap-style). Default ON. Set FALSE to revert to per-panel
    # theme-drawn axis lines (which leave a tiny visible gap between panels).
    continuous_x_axis = TRUE,
    # Knobs for the rotated CMV-status label below each panel (only active
    # when x_axis_mode = "group"). The top margin determines how far below
    # the panel bottom the label sits — needs to be large enough that the
    # rotated text's upper edge doesn't overflow back into the panel above.
    # For a 12pt-sized "Persistent IgM" at 45 deg, ~22pt is the threshold.
    x_axis_label_angle      = 45,
    # NEW (auto-anchored labels): in build_cluster_facet's group mode the
    # CMV labels are drawn as panel-scoped geom_text at the centre tick,
    # so hjust = 1 / vjust = 1 anchor the rightmost letter at the tick
    # and the rest of the label trails diagonally down-and-to-the-left.
    # (Old default was hjust = 0.5 for strip-text rendering; the new
    # geom_text approach gets the standard "tilted axis label" look only
    # with hjust = 1.)
    x_axis_label_hjust      = 1,
    x_axis_label_vjust      = 1,
    # NOTE: x_axis_label_margin_top is unused in the new geom_text path
    # (label y is computed from axis_y_value / tick_y_end). Kept for
    # backward compatibility -- it still works in inline / non-group
    # modes where the original strip-text approach renders the labels.
    x_axis_label_margin_top = 22,
    # Top-level layout. "grid" is the original 2D layout (rows = clusters,
    # cols = groups). "single_row" puts ALL clusters on one row as side-by-
    # side "blocks", each block = the cluster's 4 group panels with a
    # facet_grid-style header text above naming the cluster. Independent y
    # per cluster block; shared y within each block (when shared_y_per_row).
    # Layout values:
    #   "single_row" -- modern faceted path (uses build_cluster_facet).
    #     All cluster blocks placed in ONE row via wrap_plots(nrow = 1).
    #     Use when you have a small number of clusters (e.g. merging_act
    #     with 3 clusters fits comfortably in one row).
    #   "wrap"       -- modern faceted path. Wraps cluster blocks to
    #     multiple rows at n_cols per row via wrap_plots(ncol = n_cols).
    #     Use for many-cluster figures (e.g. merging_combined with
    #     ~15-45 clusters).
    #   "grid"       -- LEGACY per-(cluster, group) make_composite path.
    #     Doesn't go through build_cluster_facet. n_cols is IGNORED
    #     (hardcoded to length(group_levels)). Don't use for new work.
    layout            = c("grid", "single_row", "wrap"),
    # Cluster header (single_row mode only). NULL -> derived from title_size.
    cluster_header_size  = NULL,
    cluster_header_face  = "bold",
    cluster_block_gap_pt = 14,   # gap between adjacent cluster blocks
    # X-axis labelling style:
    #   "peptide" -> show peptide ticks (No Peptide / IE1 / pp65) on each
    #                panel, group_col level appears once as the strip-top
    #                title of the top row.
    #   "group"   -> hide peptide ticks; show the group_col level (e.g. CMV
    #                status) as the x-axis title below each panel. The
    #                strip-top title is also hidden (avoids double-labelling).
    #                Peptide identity is conveyed by the fill legend only.
    x_axis_mode       = c("peptide", "group"),
    # Y-axis label root
    y_axis_root       = "Proportion [%]",
    # Methods annotation (figure-legend-style title at the top of the
    # patchwork). Three modes:
    #   NULL  -> no title (silent)
    #   TRUE  -> auto-generated, e.g.
    #            "Paired Wilcoxon signed-rank vs 'No Peptide'; per-cluster BH"
    #            (test type + reference + correction + star key)
    #   "<string>" -> use that exact string verbatim
    stats_title       = TRUE,
    stats_title_size  = 11,
    stats_title_face  = "italic",
    # Diagnostic
    verbose           = TRUE
) {
    # ---- Package dependencies --------------------------------------------
    if (!requireNamespace("patchwork", quietly = TRUE))
        stop("patchwork is required: install.packages('patchwork').")
    if (!requireNamespace("ggh4x", quietly = TRUE))
        stop("ggh4x is required: install.packages('ggh4x').")
    if (!requireNamespace("rstatix", quietly = TRUE))
        stop("rstatix is required (wilcox_test): install.packages('rstatix').")
    if (isTRUE(use_richtext) && !requireNamespace("ggtext", quietly = TRUE))
        stop("ggtext is required when use_richtext = TRUE: install.packages('ggtext').")

    if (is.null(row_label_size)) row_label_size <- title_size
    x_axis_mode        <- match.arg(x_axis_mode)
    layout             <- match.arg(layout)
    p_adjust_method    <- match.arg(p_adjust_method)
    log_axis_top_style <- match.arg(log_axis_top_style)
    if (is.null(cluster_header_size)) cluster_header_size <- title_size + 1

    correction <- match.arg(correction)
    if (identical(label_col, "p.adj.signif")) {
        label_col <- switch(correction,
                            per_cluster = "p.adj.signif",
                            global      = "p.adj_global.signif")
    }

    # ---- Log-axis tick label formatter -----------------------------------
    # R's base `format()` on a numeric VECTOR pads every element to a common
    # decimal width -- so c(0.01, 0.1, 1, 10, 100) becomes "0.01", "0.10",
    # "1.00", "10.00", "100.00". For log axes we want each tick to display
    # at its natural precision: "0.01", "0.1", "1", "10", "100".
    #
    # Strategy: format each value INDIVIDUALLY (no cross-value padding),
    # then strip trailing zeros that sit after a decimal point. Crucially
    # we don't touch trailing zeros in integers like "100" -- the regex only
    # bites when there's a "." in the string.
    .fmt_log_break <- function(x) {
        vapply(x, function(v) {
            s <- format(v, scientific = FALSE, trim = TRUE)
            s <- sub("(\\..*?)0+$", "\\1", s)  # ".10" -> ".1"; ".00" -> "."
            s <- sub("\\.$",        "",   s)   # trailing "." -> ""
            s
        }, character(1))
    }

    # ---- Cluster ids -----------------------------------------------------
    # Mirror plotAbundancesA1's logic: if `k` is a colData column, use it;
    # else go through CATALYST .wl_cluster_ids().
    if (!is.null(merging_col) || k %in% names(SingleCellExperiment::colData(x))) {
        cluster_id_vec <- x[[k]]
    } else {
        cluster_id_vec <- .wl_cluster_ids(x, .wl_check_k(x, k))
    }

    # ---- Build per-donor proportion table -------------------------------
    # Per-(sample_id, cluster) cell count -> per-sample cluster proportion
    cd <- as.data.frame(SingleCellExperiment::colData(x))
    # CAPTURE original factor levels (if any) for cluster ordering.
    # When x[[k]] is a factor with explicit levels -- e.g.
    #     scefilt$merging_act <- factor(..., levels = c("Non-responder",
    #         "Non-specific", "AIM"))
    # -- those levels are the user's intended left-to-right order. We
    # save them here and use them as the default cluster order below
    # (UNLESS the caller passes an explicit `clusters_order`). Without
    # this, line 225's `as.character()` drops the factor info and the
    # downstream "unique(as.character(...))" fallback orders clusters by
    # FIRST APPEARANCE in the data, which is not what the user wants.
    .cluster_id_levels <- if (is.factor(cluster_id_vec)) {
        levels(cluster_id_vec)
    } else NULL
    cd$.cluster_id <- as.character(cluster_id_vec)

    needed_cols <- c("sample_id", patient_col, group_col, peptide_col)
    missing_cols <- setdiff(needed_cols, names(cd))
    if (length(missing_cols))
        stop("Missing required colData columns: ", paste(missing_cols, collapse = ", "))

    # CRITICAL: group_by + summarise only creates rows for (cluster, sample)
    # combinations that actually have cells. Samples with ZERO cells of a
    # particular cluster get silently DROPPED from counts_df, which means
    # they never get a freq = 0 row and never plot as a dot at miny on log
    # axes. Use tidyr::complete() (or table()-style densification) to fill
    # in the missing combinations with n_cells_cluster = 0. This matches
    # plotAbundancesA1's `table(cluster_id, sample_id)` which produces a
    # full contingency table including zero cells.

    .within_cluster_mode <- !is.null(act_col) && !is.null(aim_level)

    if (.within_cluster_mode) {
        # ---- "Within-cluster" denominator mode --------------------------
        # numerator = cells where act_col == aim_level (within each cluster
        #             panel and sample)
        # denominator = ALL cells in that cluster panel (within sample)
        if (!act_col %in% names(cd)) {
            stop("plotAimPaired (within-cluster mode): act_col '", act_col,
                 "' not found in colData.")
        }
        if (!any(cd[[act_col]] == aim_level, na.rm = TRUE)) {
            warning("plotAimPaired (within-cluster mode): no cells have ",
                    act_col, " == '", aim_level, "'. All freqs will be 0.")
        }
        counts_df <- cd %>%
            dplyr::filter(.data[[act_col]] == aim_level) %>%
            dplyr::group_by(.cluster_id, sample_id) %>%
            dplyr::summarise(n_cells_cluster = dplyr::n(), .groups = "drop") %>%
            tidyr::complete(
                .cluster_id = unique(cd$.cluster_id),
                sample_id   = unique(cd$sample_id),
                fill        = list(n_cells_cluster = 0L)
            )
        # Per-(sample, cluster) totals — denominator is "all cells in this
        # cluster panel for this sample", NOT total cells in the sample.
        totals_df <- cd %>%
            dplyr::group_by(.cluster_id, sample_id) %>%
            dplyr::summarise(n_cells_sample = dplyr::n(), .groups = "drop")
        df <- counts_df %>%
            dplyr::left_join(totals_df, by = c(".cluster_id", "sample_id")) %>%
            dplyr::mutate(
                n_cells_sample = tidyr::replace_na(n_cells_sample, 0L),
                freq = ifelse(n_cells_sample > 0,
                              n_cells_cluster / n_cells_sample * 100,
                              0)
            )
    } else {
        # ---- Default: % of parent lineage --------------------------------
        counts_df <- cd %>%
            dplyr::group_by(.cluster_id, sample_id) %>%
            dplyr::summarise(n_cells_cluster = dplyr::n(), .groups = "drop") %>%
            tidyr::complete(
                .cluster_id = unique(cd$.cluster_id),
                sample_id   = unique(cd$sample_id),
                fill        = list(n_cells_cluster = 0L)
            )

        # Total cells per sample (denominator for proportion)
        totals_df <- cd %>%
            dplyr::group_by(sample_id) %>%
            dplyr::summarise(n_cells_sample = dplyr::n(), .groups = "drop")

        df <- counts_df %>%
            dplyr::left_join(totals_df, by = "sample_id") %>%
            dplyr::mutate(freq = n_cells_cluster / n_cells_sample * 100)
    }

    # Attach metadata (patient, group, peptide) by sample
    sample_meta <- cd %>%
        dplyr::select(sample_id, dplyr::all_of(c(patient_col, group_col, peptide_col))) %>%
        dplyr::distinct()
    df <- df %>% dplyr::left_join(sample_meta, by = "sample_id")

    # ---- Filter clusters -------------------------------------------------
    if (!is.null(excluded_clusters)) {
        df <- df[!df$.cluster_id %in% excluded_clusters, , drop = FALSE]
    }
    if (!is.null(keep_clusters)) {
        # Whitelist filter. Applied AFTER proportions are calculated (above),
        # so denominators are unchanged -- this just suppresses panels.
        .missing <- setdiff(keep_clusters, unique(df$.cluster_id))
        if (length(.missing) > 0) {
            warning("keep_clusters contains cluster ids not present in data: ",
                    paste(.missing, collapse = ", "))
        }
        df <- df[df$.cluster_id %in% keep_clusters, , drop = FALSE]
    }
    if (!is.null(clusters_order)) {
        # Explicit user-supplied ordering wins.
        df$.cluster_id <- factor(df$.cluster_id, levels = clusters_order)
    } else if (!is.null(keep_clusters)) {
        # No explicit order, but keep_clusters was supplied -- use its order
        # as the panel layout (intersect drops names that weren't in data).
        df$.cluster_id <- factor(
            df$.cluster_id,
            levels = intersect(keep_clusters, unique(as.character(df$.cluster_id)))
        )
    } else if (!is.null(.cluster_id_levels)) {
        # Use the factor levels supplied via x[[k]] (the user's
        # `scefilt$merging_act <- factor(..., levels = ...)` call).
        # `intersect()` drops any unused levels while preserving order.
        df$.cluster_id <- factor(
            df$.cluster_id,
            levels = intersect(.cluster_id_levels, unique(as.character(df$.cluster_id)))
        )
    } else {
        df$.cluster_id <- factor(df$.cluster_id, levels = unique(as.character(df$.cluster_id)))
    }
    clusters <- levels(droplevels(df$.cluster_id))

    # ---- Order group_col + peptide_col -----------------------------------
    if (is.null(group_levels)) {
        group_levels <- if (is.factor(cd[[group_col]])) {
            levels(droplevels(cd[[group_col]]))
        } else {
            sort(unique(as.character(cd[[group_col]])))
        }
    }
    df[[group_col]] <- factor(df[[group_col]], levels = group_levels)

    if (is.null(peptide_levels)) {
        # Put reference first
        all_peps <- if (is.factor(cd[[peptide_col]])) {
            levels(droplevels(cd[[peptide_col]]))
        } else {
            sort(unique(as.character(cd[[peptide_col]])))
        }
        peptide_levels <- c(reference_level, setdiff(all_peps, reference_level))
    }
    df[[peptide_col]] <- factor(df[[peptide_col]], levels = peptide_levels)

    if (!(reference_level %in% peptide_levels))
        stop("reference_level '", reference_level, "' not found in peptide_col levels.")

    # ---- Paired Wilcoxon stats per (cluster, group, non-reference peptide) ----
    non_ref_peps <- setdiff(peptide_levels, reference_level)
    full_stats <- NULL
    if (isTRUE(show_stats) && length(non_ref_peps) >= 1) {
        # For each (cluster, group), test reference vs each non-reference peptide
        # via paired Wilcoxon. Pairing is by patient_col; donors missing either
        # the reference or the peptide sample are dropped from that single test.
        stat_rows <- list()
        for (cl in clusters) for (gl in group_levels) for (pep in non_ref_peps) {
            sub <- df[df$.cluster_id == cl &
                      as.character(df[[group_col]]) == gl &
                      as.character(df[[peptide_col]]) %in% c(reference_level, pep),
                      , drop = FALSE]
            # Pivot to wide so wilcox_test can match pairs by patient_col
            wide <- sub %>%
                dplyr::select(dplyr::all_of(c(patient_col, peptide_col)), freq) %>%
                tidyr::pivot_wider(names_from = !!rlang::sym(peptide_col),
                                   values_from = freq)
            # Guard: if all rows for one peptide were absent from `sub`,
            # pivot_wider won't create that column. wide[[<missing>]] returns
            # NULL and is.na(NULL) is logical(0), which can't be used as a
            # row mask. In that case there are zero usable pairs -- emit an
            # NA stat row and skip the test.
            if (!(reference_level %in% names(wide)) || !(pep %in% names(wide))) {
                stat_rows[[length(stat_rows) + 1L]] <- data.frame(
                    .cluster_id = cl, group = gl,
                    group1 = reference_level, group2 = pep,
                    n1 = 0L, n2 = 0L,
                    statistic = NA_real_, p = NA_real_,
                    stringsAsFactors = FALSE
                )
                next
            }
            wide <- wide[!is.na(wide[[reference_level]]) & !is.na(wide[[pep]]), , drop = FALSE]
            n_pairs <- nrow(wide)
            if (n_pairs < 2) {
                stat_rows[[length(stat_rows) + 1L]] <- data.frame(
                    .cluster_id = cl, group = gl, group1 = reference_level, group2 = pep,
                    n1 = n_pairs, n2 = n_pairs,
                    statistic = NA_real_, p = NA_real_,
                    stringsAsFactors = FALSE
                )
                next
            }
            wt <- tryCatch(
                stats::wilcox.test(wide[[reference_level]], wide[[pep]], paired = TRUE,
                                   exact = FALSE),
                error = function(e) NULL
            )
            stat_rows[[length(stat_rows) + 1L]] <- data.frame(
                .cluster_id = cl, group = gl,
                group1 = reference_level, group2 = pep,
                n1 = n_pairs, n2 = n_pairs,
                statistic = if (is.null(wt)) NA_real_ else unname(wt$statistic),
                p = if (is.null(wt)) NA_real_ else wt$p.value,
                stringsAsFactors = FALSE
            )
        }
        full_stats <- do.call(rbind, stat_rows)

        # Per-cluster correction (over each cluster's group_levels x non_ref_peps
        # tests). In this function there's no "per-pair" correction layer (every
        # test is paired against the same reference within one cluster x group
        # cell), so p.adj IS the per-cluster correction. We expose it under both
        # names (`p.adj` and `p.adj_per_cluster`) so callers can use the strip-
        # plotter convention `label_col = "p.adj_per_cluster.signif"` and it
        # Just Works. Method = `p_adjust_method` (default "holm").
        full_stats <- full_stats %>%
            dplyr::group_by(.cluster_id) %>%
            dplyr::mutate(p.adj = stats::p.adjust(p, method = p_adjust_method)) %>%
            dplyr::ungroup() %>%
            as.data.frame()
        full_stats$p.adj_per_cluster <- full_stats$p.adj   # alias

        # Global correction (over the whole figure). Same method as per-cluster
        # for consistency in the figure-legend description.
        full_stats$p.adj_global <- stats::p.adjust(full_stats$p, method = p_adjust_method)

        # Signif columns
        .signif <- function(v) ifelse(is.na(v), "ns",
                                      dplyr::case_when(
                                          v < 0.0001 ~ "****",
                                          v < 0.001  ~ "***",
                                          v < 0.01   ~ "**",
                                          v < 0.05   ~ "*",
                                          TRUE       ~ "ns"))
        full_stats$p.adj.signif             <- .signif(full_stats$p.adj)
        full_stats$p.adj_per_cluster.signif <- full_stats$p.adj.signif  # alias
        full_stats$p.adj_global.signif      <- .signif(full_stats$p.adj_global)
        full_stats$signif_changed <- full_stats$p.adj.signif != full_stats$p.adj_global.signif
        # Reorder columns (kept verbose -- mirrors the other strip plotters)
        full_stats <- full_stats[, c(".cluster_id", "group", "group1", "group2",
                                      "n1", "n2", "statistic", "p",
                                      "p.adj", "p.adj.signif",
                                      "p.adj_per_cluster", "p.adj_per_cluster.signif",
                                      "p.adj_global", "p.adj_global.signif",
                                      "signif_changed")]
    }

    # ---- Validate label_col -------------------------------------------------
    # If the caller asked for a column we don't produce (e.g. a typo or a
    # column from an external_stats convention this function doesn't carry),
    # fall back to p.adj.signif with a warning rather than failing inside
    # build_strip_plot with a confusing data.frame length-mismatch error.
    if (isTRUE(show_stats) && !is.null(full_stats)) {
        if (!(label_col %in% names(full_stats))) {
            warning("plotAimPaired: label_col '", label_col,
                    "' not in stats columns (",
                    paste(names(full_stats), collapse = ", "),
                    "). Falling back to 'p.adj.signif'.", call. = FALSE)
            label_col <- "p.adj.signif"
        }
    }

    # ---- Label transform (asterisk styling, same convention as siblings) ----
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

    # ---- Strip / data height fractions (uniform across cells) ------------
    strip_h <- strip_height_fraction
    data_h  <- 1 - strip_height_fraction
    strip_aspect <- facet_ratio * strip_h / data_h

    # ---- Per-cell stats lookup (cluster x group key) ---------------------
    get_cell_stats <- function(cl, gl) {
        if (is.null(full_stats)) return(NULL)
        sub <- full_stats[full_stats$.cluster_id == cl & full_stats$group == gl, , drop = FALSE]
        if (isTRUE(hide_ns)) {
            sub <- sub[as.character(sub[[label_col]]) != "ns", , drop = FALSE]
        }
        sub
    }
    max_brackets <- length(non_ref_peps)
    if (isTRUE(hide_ns) && !is.null(full_stats)) {
        # cap max_brackets at the worst-case (most surviving in any single cell)
        cap <- 0L
        for (cl in clusters) for (gl in group_levels) {
            sub <- full_stats[full_stats$.cluster_id == cl & full_stats$group == gl, , drop = FALSE]
            sub <- sub[as.character(sub[[label_col]]) != "ns", , drop = FALSE]
            cap <- max(cap, nrow(sub))
        }
        max_brackets <- cap
    }

    # ---- Builders --------------------------------------------------------
    build_data_plot <- function(df_cell, cl, gl,
                                show_row_label, show_col_title,
                                force_y_top = NULL,
                                include_cluster_in_y_label = TRUE) {
        # In `single_row` layout the cluster name lives in the header above
        # the block (facet_grid-style), NOT in the y-axis title. So we
        # optionally suppress the cluster prefix on the y label.
        y_lab <- if (isTRUE(show_row_label)) {
            if (isTRUE(include_cluster_in_y_label)) paste0(cl, "\n", y_axis_root) else y_axis_root
        } else NULL

        # Pre-shift for log so zeros become miny
        df_cell$y_plot <- if (isTRUE(log)) df_cell$freq + miny else df_cell$freq

        # ---- Pre-compute x positions (peptide factor position + jitter) ---
        # Single source of truth for x. Lines, points, and downstream layers
        # all read .x_pos from df_cell. This is the same trick used in
        # plotAbundancesA1's connect_donors path: precomputing avoids the
        # ggbeeswarm / position_dodge invisible-position trap that mis-aligns
        # line endpoints with points.
        pep_factor <- factor(as.character(df_cell[[peptide_col]]),
                             levels = peptide_levels)
        df_cell$.x_base <- as.numeric(pep_factor)
        if (jitter_width > 0 && nrow(df_cell) > 0) {
            old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
                           get(".Random.seed", envir = .GlobalEnv) else NULL
            set.seed(jitter_seed)
            df_cell$.x_jit <- runif(nrow(df_cell), -jitter_width, jitter_width)
            if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
        } else {
            df_cell$.x_jit <- 0
        }
        df_cell$.x_pos <- df_cell$.x_base + df_cell$.x_jit

        # ---- Build fan-from-reference donor lines -----------------------
        # IE1 and pp65 are PARALLEL stim conditions, not sequential. So per
        # donor we draw separate lines from the reference (No Peptide) to
        # each non-reference peptide -- NOT a single polyline through all
        # three points. A donor missing either endpoint of a given pair is
        # dropped from THAT pair only. fan_df rows carry .x_pos from df_cell
        # so endpoints land exactly on the points.
        .build_fan_lines <- function(df_in) {
            non_ref <- setdiff(peptide_levels, reference_level)
            out <- list()
            for (pep in non_ref) {
                sub <- df_in[as.character(df_in[[peptide_col]])
                             %in% c(reference_level, pep), , drop = FALSE]
                # Need exactly 2 points (the reference and `pep`) per donor.
                keep_donors <- names(table(sub[[patient_col]]))[
                    table(sub[[patient_col]]) == 2
                ]
                sub <- sub[as.character(sub[[patient_col]]) %in% keep_donors, , drop = FALSE]
                if (nrow(sub) > 0) {
                    sub$.pair_id <- paste0(as.character(sub[[patient_col]]), "__", pep)
                    out[[pep]] <- sub
                }
            }
            if (length(out) == 0) return(NULL)
            do.call(rbind, out)
        }
        fan_df <- .build_fan_lines(df_cell)

        # Base ggplot. fill = peptide_col is inherited; geom_boxplot uses
        # it for column fills, geom_point uses it for fill of the dots.
        # Lines override with their own grey color (no inheritance).
        p <- ggplot2::ggplot(
                df_cell,
                ggplot2::aes(x    = .data[[peptide_col]],
                             y    = y_plot,
                             fill = .data[[peptide_col]])
            )

        # Boxplot underlay (toggleable). Uses the discrete peptide_col x so
        # the box sits at the column center; alpha < 1 so points / lines
        # passing through the IQR remain visible.
        if (isTRUE(show_boxplots)) {
            p <- p + ggplot2::geom_boxplot(
                color         = "black",
                alpha         = boxplot_alpha,
                linewidth     = boxplot_linewidth,
                width         = boxplot_width,
                outlier.color = NA,
                show.legend   = FALSE
            )
        }

        # Donor fan-lines (drawn behind points). Each .pair_id is a 2-point
        # group: (reference, one non-reference peptide). x = .x_pos so the
        # line connects the SAME jittered x values the points sit at.
        if (!is.null(fan_df)) {
            p <- p + ggplot2::geom_line(
                data    = fan_df,
                mapping = ggplot2::aes(x = .x_pos, y = y_plot,
                                       group = .pair_id),
                color = line_color, linewidth = line_size, alpha = line_alpha,
                inherit.aes = FALSE,
                show.legend = FALSE
            )
        }
        # Points: explicit x = .x_pos (no geom_quasirandom, no dodge) so the
        # line endpoints land exactly on each dot.
        p <- p +
            ggplot2::geom_point(
                mapping = ggplot2::aes(x = .x_pos, y = y_plot,
                                       fill = .data[[peptide_col]]),
                inherit.aes = FALSE,
                shape = 21, color = "black", stroke = 0.4, size = point_size,
                show.legend = TRUE
            ) +
            ggplot2::labs(
                # x label depends on x_axis_mode (set in the assembly below):
                # "peptide" -> NULL (peptide names are the per-tick text instead)
                # "group"   -> the group_col level for this panel (e.g. "Never")
                x = if (identical(x_axis_mode, "group")) gl else NULL,
                y = y_lab,
                title = if (isTRUE(show_col_title) && !isTRUE(title_above_strip)) gl else NULL
            ) +
            ggplot2::scale_x_discrete(limits = peptide_levels, drop = FALSE) +
            ggh4x::force_panelsizes(cols = grid::unit(panel_width_cm, "cm")) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                panel.grid       = ggplot2::element_blank(),
                aspect.ratio     = facet_ratio,
                plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5,
                                                         size = title_size,
                                                         margin = ggplot2::margin(b = 2)),
                axis.title.y     = ggplot2::element_text(size = row_label_size, face = "bold"),
                # x-axis title is used in `group` mode to carry the CMV-status
                # label below each panel. All four positioning knobs are
                # user-tunable (x_axis_label_angle / hjust / vjust /
                # margin_top). With angle = 45, hjust = 0.5, vjust = 1, the
                # rotated text's top-center anchors at the panel's centered
                # tick position. The top margin determines how far below the
                # panel the label sits — needs to be ≳ half_text_width *
                # sin(angle) so the rotated upper-right doesn't overflow.
                axis.title.x     = ggplot2::element_text(
                                       size   = axis_text_size,
                                       face   = "plain",
                                       angle  = x_axis_label_angle,
                                       hjust  = x_axis_label_hjust,
                                       vjust  = x_axis_label_vjust,
                                       margin = ggplot2::margin(t = x_axis_label_margin_top,
                                                                b = 2)
                                   ),
                axis.text        = ggplot2::element_text(color = "black", size = axis_text_size),
                axis.text.x      = if (identical(x_axis_mode, "group")) {
                                       ggplot2::element_blank()
                                   } else {
                                       ggplot2::element_text(angle = 45, hjust = 1, vjust = 1,
                                                             size = axis_text_size)
                                   },
                axis.ticks.x     = if (identical(x_axis_mode, "group"))
                                       ggplot2::element_blank()
                                   else ggplot2::element_line(color = "black"),
                axis.ticks       = ggplot2::element_line(color = "black"),
                plot.margin      = ggplot2::margin(2, 2, 2, 2)
            )

        # ---- Y-axis scaling
        # `force_y_top` (when non-NULL) is the pre-computed row-shared max
        # for this cluster (data max across every group in the row). Used
        # only when shared_y_per_row = TRUE in the assembly below.
        if (isTRUE(log)) {
            decade_breaks <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)
            # "Nice" log breaks at log-2-5 spacing — used as termination
            # candidates so the axis doesn't waste a whole decade above sparse
            # data (y_max = 1.2 -> axis_top = 2, not 10).
            nice_log_breaks <- c(
                1e-4, 2e-4, 5e-4,
                1e-3, 2e-3, 5e-3,
                1e-2, 2e-2, 5e-2,
                1e-1, 2e-1, 5e-1,
                1,    2,    5,
                10,   20,   50,
                100,  200,  500,
                1000, 2000, 5000, 10000
            )

            y_max <- if (!is.null(force_y_top)) force_y_top
                     else                       max(df_cell$y_plot, na.rm = TRUE)
            if (!is.finite(y_max) || y_max <= 0) y_max <- 1

            target <- y_max * (1 + tight_top_axis_overhang)

            if (identical(log_axis_top_style, "decade")) {
                # Legacy behaviour: snap to next power of 10
                axis_top <- 10 ^ ceiling(log10(target))
            } else {
                # Smallest nice_log_break above target
                idx <- which(nice_log_breaks >= target)[1]
                axis_top <- if (length(idx) && !is.na(idx)) {
                    nice_log_breaks[idx]
                } else {
                    10 ^ ceiling(log10(target))   # fallback if data exceeds the table
                }
            }

            # Axis ticks: decade ticks within [miny, axis_top], plus axis_top
            # itself if it's an intermediate (2x or 5x) value.
            breaks_vec <- decade_breaks[decade_breaks >= miny & decade_breaks <= axis_top]
            if (!(axis_top %in% breaks_vec)) {
                breaks_vec <- sort(c(breaks_vec, axis_top))
            }

            p <- p + ggplot2::scale_y_continuous(
                trans  = "log10",
                limits = c(miny, axis_top),
                breaks = breaks_vec,
                labels = .fmt_log_break(breaks_vec)
            ) +
            # clip = "off" so points sitting exactly at the bottom (`miny`,
            # which is where zeros get shifted to via the pseudocount) still
            # render fully instead of being clipped by the panel boundary.
            # Matches plotAbundancesA1's log path.
            ggplot2::coord_cartesian(clip = "off")
        } else if (isTRUE(tight_top_axis)) {
            y_max <- if (!is.null(force_y_top)) force_y_top
                     else                       max(df_cell$y_plot, na.rm = TRUE)
            if (!is.finite(y_max) || y_max <= 0) y_max <- 1
            breaks_vec <- pretty(c(0, y_max * 1.001), n = 5)
            top_break  <- max(breaks_vec)
            axis_top   <- top_break * (1 + tight_top_axis_overhang)
            p <- p + ggplot2::scale_y_continuous(
                breaks = breaks_vec,
                limits = c(0, axis_top),
                expand = ggplot2::expansion(mult = c(y_expand_low_mult, 0))
            )
        } else {
            p <- p + ggplot2::scale_y_continuous(
                limits = if (!is.null(force_y_top)) c(0, force_y_top)
                         else                       c(0, NA),
                expand = ggplot2::expansion(mult = c(y_expand_low_mult, y_expand_high_mult))
            )
        }

        # ---- Hide y-axis chrome on non-leftmost panels (when shared_y_per_row).
        # `show_row_label` is the leftmost-column flag; suppressing y axis
        # text/ticks/line on every other panel matches the facet_wrap look
        # (axis appears once on the left of each row).
        if (!isTRUE(show_row_label)) {
            p <- p + ggplot2::theme(
                axis.text.y  = ggplot2::element_blank(),
                axis.ticks.y = ggplot2::element_blank(),
                axis.line.y  = ggplot2::element_blank(),
                axis.title.y = ggplot2::element_blank()
            )
        }

        if (!is.null(fill_palette)) {
            p <- p + ggplot2::scale_fill_manual(values = fill_palette, name = peptide_col, drop = FALSE)
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

        # NOTE: `continuous_x_axis` parameter is currently a no-op while we
        # work out a clean cross-panel axis approach. Earlier attempts:
        #   - annotate("segment", x = numeric, ...): corrupted the discrete
        #     x scale by trying to add the numeric values as new levels.
        #   - annotate("segment", x = -Inf/Inf, y = -Inf/Inf, ...): worked
        #     on linear scales but on log scales the y = -Inf endpoint was
        #     pushed through log10() -> NaN -> the segment was dropped as
        #     "outside scale range" with warnings.
        # The proper fix likely requires either (a) a refactor to use
        # facet_grid2/facet_wrap2 within each cluster block (so a single
        # shared axis line spans the block), or (b) a grid::grob-based
        # cross-panel annotation via patchwork::inset_element. Leaving the
        # theme's axis.line.x in place for now so at least each panel has
        # its own clean axis line.

        p
    }

    build_strip_plot <- function(stat_cell, cl, gl, show_col_title) {
        n_brackets <- if (is.null(stat_cell)) 0 else nrow(stat_cell)
        n_x        <- length(peptide_levels)
        y_axis_max <- max(max_brackets, 1)

        if (n_brackets > 0) {
            bracket_y <- n_brackets - seq_len(n_brackets) + 0.5
            x1 <- match(as.character(stat_cell$group1), peptide_levels)
            x2 <- match(as.character(stat_cell$group2), peptide_levels)
            x_mid <- (x1 + x2) / 2
            raw_labels <- as.character(stat_cell[[label_col]])
            labels     <- label_transform(raw_labels)

            seg_h <- data.frame(x = x1, xend = x2, y = bracket_y, yend = bracket_y)
            seg_l <- data.frame(x = x1, xend = x1, y = bracket_y, yend = bracket_y - bracket_tip_npc)
            seg_r <- data.frame(x = x2, xend = x2, y = bracket_y, yend = bracket_y - bracket_tip_npc)
            seg_df <- rbind(seg_h, seg_l, seg_r)

            is_ast <- grepl("[*]", raw_labels)
            label_y <- bracket_y + label_nudge_npc
            label_y[is_ast] <- label_y[is_ast] + asterisk_y_offset
            vj_vec <- ifelse(is_ast, asterisk_vjust, label_vjust)
            text_df <- data.frame(x = x_mid, y = label_y, label = labels, vj = vj_vec)
        } else {
            seg_df  <- data.frame(x = numeric(0), xend = numeric(0),
                                  y = numeric(0), yend = numeric(0))
            text_df <- data.frame(x = numeric(0), y = numeric(0),
                                  label = character(0), vj = numeric(0))
        }

        p <- ggplot2::ggplot()
        if (nrow(seg_df) > 0) {
            p <- p + ggplot2::geom_segment(
                data    = seg_df,
                mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                color   = bracket_color, linewidth = bracket_size,
                lineend = "round", linejoin = "round"
            )
            if (isTRUE(use_richtext)) {
                p <- p + ggtext::geom_richtext(
                    data    = text_df,
                    mapping = ggplot2::aes(x = x, y = y, label = label, vjust = vj),
                    size    = label_size, fontface = label_fontface,
                    label.size = NA, fill = NA,
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
            # Hide strip-top title in `group` x_axis_mode to avoid double-
            # labelling (the same CMV status would otherwise appear both
            # above the strip AND below the panel).
            ggplot2::labs(title = if (isTRUE(show_col_title) && isTRUE(title_above_strip) &&
                                       !identical(x_axis_mode, "group")) gl else NULL) +
            ggplot2::scale_x_continuous(limits = c(0.5, n_x + 0.5), expand = c(0, 0)) +
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
        p
    }

    # ---- Legend appearance helpers --------------------------------------
    # Applied via `&` to the outer patchwork so they propagate to all
    # inner ggplots' legends, which patchwork's `guides = "collect"` then
    # merges into a single right-side legend.
    .build_legend_theme <- function(pos, key_sz_mm, txt_sz, ttl_sz, vsp_pt,
                                     ax_txt, row_lbl) {
        if (identical(pos, "none")) return(ggplot2::theme(legend.position = "none"))
        eff_txt <- if (is.null(txt_sz)) ax_txt   else txt_sz
        eff_ttl <- if (is.null(ttl_sz)) row_lbl  else ttl_sz
        ggplot2::theme(
            legend.position     = pos,
            legend.key.size     = grid::unit(key_sz_mm, "mm"),
            legend.key.height   = grid::unit(key_sz_mm, "mm"),
            legend.key.width    = grid::unit(key_sz_mm, "mm"),
            legend.text         = ggplot2::element_text(size = eff_txt),
            legend.title        = ggplot2::element_text(size = eff_ttl, face = "bold"),
            # Tighten vertical whitespace BETWEEN legend entries.
            legend.spacing.y    = grid::unit(vsp_pt, "pt"),
            # Required for legend.spacing.y to actually apply (ggplot2 quirk):
            # the legend layout collapses spacing.y unless byrow is forced
            # via guide_legend()'s nrow / byrow; we set byrow via the guide
            # helper below. Here we also tighten the per-key spacing.
            legend.key          = ggplot2::element_blank()
        )
    }
    # `override.aes = list(size = ...)` forces the legend GLYPH (the
    # filled circle) to render at the new larger size, independent of the
    # actual point_size used inside each panel.
    .build_fill_guide <- function(fill_var_name, key_sz_mm) {
        # Glyph "size" is in mm-pt geom units. Translate key.size mm to
        # geom size: ~ key_sz_mm * 1.4 for visually-matched circle area.
        glyph_size <- key_sz_mm * 1.4
        ggplot2::guides(
            fill = ggplot2::guide_legend(
                override.aes = list(
                    size   = glyph_size,
                    shape  = 21,
                    stroke = 0.6,
                    color  = "black"
                ),
                byrow = TRUE   # makes legend.spacing.y actually apply
            )
        )
    }

    # ---- Build the (cluster x group) grid of composites ------------------
    # Rows = cluster, cols = group level. Composite per cell = [strip / data]
    composites <- list()
    idx <- 0L
    n_cols_grid <- length(group_levels)

    # ---- Precompute per-cluster y-top for shared_y_per_row -------------
    # When shared_y_per_row = TRUE, every panel in one cluster's row uses
    # the same y-axis termination, computed from the cluster's max across
    # all groups. Lets us hide the y-axis chrome on non-leftmost panels
    # without producing visually wrong (misaligned) panels.
    y_top_per_cluster <- list()
    if (isTRUE(shared_y_per_row)) {
        for (cl in clusters) {
            sub_cl <- df[as.character(df$.cluster_id) == cl, , drop = FALSE]
            if (nrow(sub_cl) == 0) {
                y_top_per_cluster[[cl]] <- 1
            } else {
                y_max <- if (isTRUE(log)) {
                    max(sub_cl$freq + miny, na.rm = TRUE)
                } else {
                    max(sub_cl$freq, na.rm = TRUE)
                }
                if (!is.finite(y_max) || y_max <= 0) y_max <- 1
                y_top_per_cluster[[cl]] <- y_max
            }
        }
    }

    # ---- Helper: build one (cluster, group) composite -----------------------
    # `col_pos` is the 1-based column position WITHIN whatever row this cell
    # belongs to (for "grid" that's the global group-column; for
    # "single_row" that's the column within this cluster's 4-panel block).
    # `include_cluster_in_y_label` is FALSE in single_row mode because the
    # cluster name lives in the header above the block, not on the y-axis.
    make_composite <- function(cl, gl, col_pos,
                               include_cluster_in_y_label,
                               show_col_title_for_strip) {
        df_cell <- df[as.character(df$.cluster_id) == cl &
                      as.character(df[[group_col]]) == gl, , drop = FALSE]
        if (nrow(df_cell) == 0) {
            df_cell <- data.frame(
                .cluster_id = cl,
                freq        = c(0, 1),
                y_plot      = c(0, 1)
            )
            df_cell[[group_col]]   <- factor(gl, levels = group_levels)
            df_cell[[peptide_col]] <- factor(peptide_levels[c(1, length(peptide_levels))],
                                             levels = peptide_levels)
            df_cell[[patient_col]] <- "dummy"
        }
        stat_cell <- get_cell_stats(cl, gl)

        show_row_label <- (col_pos == 1L)
        if (isTRUE(shared_y_per_row)) {
            force_y_top <- y_top_per_cluster[[cl]]
        } else {
            force_y_top <- NULL
            show_row_label <- TRUE
        }

        p_data <- build_data_plot(
            df_cell, cl, gl,
            show_row_label = show_row_label,
            show_col_title = show_col_title_for_strip,
            force_y_top    = force_y_top,
            include_cluster_in_y_label = include_cluster_in_y_label
        )
        p_strip <- build_strip_plot(
            stat_cell, cl, gl,
            show_col_title = show_col_title_for_strip
        )

        # Per-composite plot.margin: tight LR within a row, wider TB between rows
        if (isTRUE(tighten_within_row)) {
            gap_v <- panel_gap_pt
            gap_x <- panel_gap_x_pt
            comp_margin <- if (col_pos == 1L) {
                ggplot2::margin(t = gap_v, r = gap_x, b = gap_v, l = gap_v, unit = "pt")
            } else {
                ggplot2::margin(t = gap_v, r = gap_x, b = gap_v, l = gap_x, unit = "pt")
            }
        } else {
            comp_margin <- ggplot2::margin(
                t = panel_gap_pt, r = panel_gap_pt,
                b = panel_gap_pt, l = panel_gap_pt, unit = "pt"
            )
        }

        p_strip / p_data +
            patchwork::plot_layout(heights = c(strip_h, data_h)) &
            ggplot2::theme(plot.margin = comp_margin)
    }

    # ---- Branch on layout ---------------------------------------------------
    # `single_row` and `wrap` share the same modern faceted construction
    # (build_cluster_facet per cluster), differing only in how the final
    # cluster blocks are wrapped into the figure.
    if (identical(layout, "single_row") || identical(layout, "wrap")) {
        # Each cluster is built as a single faceted ggplot via
        # ggh4x::facet_wrap2(vars(group_col), nrow = 1), with the stats
        # brackets in a sibling faceted ggplot above. Patchwork stacks them
        # [strip / data] per cluster, then concatenates cluster blocks.
        #
        # The win vs the previous per-cell composites: within one ggplot the
        # x-axis line is rendered ONCE across all facets, giving a genuinely
        # continuous axis at panel.spacing = 0. No sub-pixel patchwork gap.

        # ---- Helper: build one cluster's faceted data + stats composite ----
        build_cluster_facet <- function(cl) {
            # Cluster-level data (all groups for this cluster)
            df_cl <- df[as.character(df$.cluster_id) == cl, , drop = FALSE]
            if (nrow(df_cl) == 0) {
                # Dummy data so the facet still renders all expected panels
                df_cl <- data.frame(
                    .cluster_id = cl,
                    freq        = rep(0, length(group_levels) * length(peptide_levels))
                )
                df_cl[[group_col]]   <- factor(rep(group_levels, each = length(peptide_levels)),
                                               levels = group_levels)
                df_cl[[peptide_col]] <- factor(rep(peptide_levels, times = length(group_levels)),
                                               levels = peptide_levels)
                df_cl[[patient_col]] <- "dummy"
            }
            # Ensure factor ordering
            df_cl[[group_col]]   <- factor(as.character(df_cl[[group_col]]), levels = group_levels)
            df_cl[[peptide_col]] <- factor(as.character(df_cl[[peptide_col]]), levels = peptide_levels)

            # ---- Adaptive miny per cluster -------------------------------
            # If `auto_miny` is on and the cluster's minimum non-zero
            # frequency is >= 0.1, raise miny to 0.1 so the y axis doesn't
            # extend uselessly down to 0.01 for clusters where no data
            # comes anywhere near there. Local assignment shadows the
            # outer-scope miny within this function call (R lexical
            # scoping), so every downstream reference inside
            # build_cluster_facet -- y-shift, scale_y limits, log
            # subticks, the manual axis-line y position in group mode --
            # picks up the new value automatically.
            if (isTRUE(auto_miny) && isTRUE(log)) {
                .cluster_min_nonzero <- suppressWarnings(
                    min(df_cl$freq[df_cl$freq > 0], na.rm = TRUE)
                )
                if (is.finite(.cluster_min_nonzero) &&
                    .cluster_min_nonzero >= 0.1) {
                    miny <- 0.1
                }
            }

            # y_plot (log shift)
            df_cl$y_plot <- if (isTRUE(log)) df_cl$freq + miny else df_cl$freq

            # .x_pos = peptide factor position + deterministic jitter
            pep_fac      <- factor(as.character(df_cl[[peptide_col]]), levels = peptide_levels)
            df_cl$.x_base <- as.numeric(pep_fac)
            if (jitter_width > 0 && nrow(df_cl) > 0) {
                old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
                               get(".Random.seed", envir = .GlobalEnv) else NULL
                set.seed(jitter_seed)
                df_cl$.x_jit <- runif(nrow(df_cl), -jitter_width, jitter_width)
                if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
            } else {
                df_cl$.x_jit <- 0
            }
            df_cl$.x_pos <- df_cl$.x_base + df_cl$.x_jit

            # ---- Fan lines across all groups for this cluster --------------
            # For each non-reference peptide, pair-by-(patient, group): both
            # the reference and peptide rows must exist for the same donor
            # within the same group.
            non_ref <- setdiff(peptide_levels, reference_level)
            fan_parts <- list()
            for (pep in non_ref) {
                sub <- df_cl[as.character(df_cl[[peptide_col]]) %in% c(reference_level, pep), , drop = FALSE]
                gkey <- paste0(as.character(sub[[group_col]]),   "\037",
                               as.character(sub[[patient_col]]))
                keep <- gkey %in% names(table(gkey))[table(gkey) == 2]
                sub  <- sub[keep, , drop = FALSE]
                if (nrow(sub) > 0) {
                    sub$.pair_id <- paste0(
                        as.character(sub[[group_col]]),   "\037",
                        as.character(sub[[patient_col]]), "\037",
                        pep
                    )
                    fan_parts[[pep]] <- sub
                }
            }
            fan_df <- if (length(fan_parts)) do.call(rbind, fan_parts) else NULL

            # ---- Stats: combine all groups' brackets for this cluster ------
            cluster_stats <- if (!is.null(full_stats)) {
                full_stats[full_stats$.cluster_id == cl, , drop = FALSE]
            } else NULL
            if (!is.null(cluster_stats) && isTRUE(hide_ns)) {
                cluster_stats <- cluster_stats[
                    as.character(cluster_stats[[label_col]]) != "ns", , drop = FALSE]
            }
            # Worst-case brackets-per-group for strip height (consistent
            # vertical spacing across facets within this cluster).
            max_brackets_cl <- 0L
            if (!is.null(cluster_stats) && nrow(cluster_stats) > 0) {
                max_brackets_cl <- max(0L, table(cluster_stats$group),
                                        na.rm = TRUE)
            }
            y_axis_max_strip <- max(max_brackets_cl, 1)

            # Build seg/text per (cluster, group) and rbind with a group tag
            seg_list  <- list()
            text_list <- list()
            if (!is.null(cluster_stats) && nrow(cluster_stats) > 0) {
                for (gl in group_levels) {
                    s <- cluster_stats[cluster_stats$group == gl, , drop = FALSE]
                    n_brackets <- nrow(s)
                    if (n_brackets == 0) next
                    bracket_y <- n_brackets - seq_len(n_brackets) + 0.5
                    x1 <- match(as.character(s$group1), peptide_levels)
                    x2 <- match(as.character(s$group2), peptide_levels)
                    x_mid <- (x1 + x2) / 2
                    raw_labels <- as.character(s[[label_col]])
                    labels     <- label_transform(raw_labels)

                    seg_h <- data.frame(x = x1,   xend = x2, y = bracket_y, yend = bracket_y)
                    seg_l <- data.frame(x = x1,   xend = x1, y = bracket_y, yend = bracket_y - bracket_tip_npc)
                    seg_r <- data.frame(x = x2,   xend = x2, y = bracket_y, yend = bracket_y - bracket_tip_npc)
                    seg_df <- rbind(seg_h, seg_l, seg_r)
                    seg_df[[group_col]] <- factor(gl, levels = group_levels)
                    seg_list[[gl]] <- seg_df

                    is_ast <- grepl("[*]", raw_labels)
                    label_y <- bracket_y + label_nudge_npc
                    label_y[is_ast] <- label_y[is_ast] + asterisk_y_offset
                    vj_vec <- ifelse(is_ast, asterisk_vjust, label_vjust)
                    text_df <- data.frame(
                        x = x_mid, y = label_y, label = labels, vj = vj_vec
                    )
                    text_df[[group_col]] <- factor(gl, levels = group_levels)
                    text_list[[gl]] <- text_df
                }
            }
            seg_df_all  <- if (length(seg_list))  do.call(rbind, seg_list)  else NULL
            text_df_all <- if (length(text_list)) do.call(rbind, text_list) else NULL

            # ---- Data ggplot with facet_wrap2 ------------------------------
            # y-axis range: shared across the cluster's facets via fixed
            # scale_y. tight_top axis behaviour reuses our log-2-5 / pretty
            # tick selection on the cluster-level max.
            force_y_top <- if (isTRUE(shared_y_per_row)) y_top_per_cluster[[cl]]
                            else                          max(df_cl$y_plot, na.rm = TRUE)
            if (!is.finite(force_y_top) || force_y_top <= 0) force_y_top <- 1

            p_data <- ggplot2::ggplot(
                df_cl,
                ggplot2::aes(x    = .data[[peptide_col]],
                             y    = y_plot,
                             fill = .data[[peptide_col]])
            )
            if (isTRUE(show_boxplots)) {
                p_data <- p_data + ggplot2::geom_boxplot(
                    color = "black", alpha = boxplot_alpha,
                    linewidth = boxplot_linewidth, width = boxplot_width,
                    outlier.color = NA, show.legend = FALSE
                )
            }
            if (!is.null(fan_df) && nrow(fan_df) > 0) {
                p_data <- p_data + ggplot2::geom_line(
                    data    = fan_df,
                    mapping = ggplot2::aes(x = .x_pos, y = y_plot,
                                           group = .pair_id),
                    color = line_color, linewidth = line_size, alpha = line_alpha,
                    inherit.aes = FALSE, show.legend = FALSE
                )
            }
            p_data <- p_data + ggplot2::geom_point(
                mapping = ggplot2::aes(x = .x_pos, y = y_plot,
                                       fill = .data[[peptide_col]]),
                inherit.aes = FALSE,
                shape = 21, color = "black", stroke = 0.4,
                size = point_size, show.legend = TRUE
            )

            # Y-scale (log vs linear), with force_y_top
            if (isTRUE(log)) {
                decade_breaks <- c(0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000)
                nice_log_breaks <- c(
                    1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3,
                    1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1,
                    1, 2, 5, 10, 20, 50, 100, 200, 500,
                    1000, 2000, 5000, 10000
                )
                target <- force_y_top * (1 + tight_top_axis_overhang)
                axis_top <- if (identical(log_axis_top_style, "decade")) {
                    10 ^ ceiling(log10(target))
                } else {
                    idx <- which(nice_log_breaks >= target)[1]
                    if (length(idx) && !is.na(idx)) nice_log_breaks[idx]
                    else                              10 ^ ceiling(log10(target))
                }
                # Only DECADE ticks on the axis — don't add the intermediate
                # axis_top (e.g. 200) as a tick label, which previously made
                # the y-axis sprout odd-looking "200" / "500" / "20" labels.
                # The axis line still extends to axis_top, just with no tick
                # there.
                breaks_vec <- decade_breaks[decade_breaks >= miny & decade_breaks <= axis_top]

                # ---- All breaks (decades + subticks) with empty labels
                # for subticks. Using ONE set of breaks via scale_y means
                # the theme's axis.ticks.y mechanism draws them, which
                # automatically respects facet behavior (renders only on
                # the leftmost facet in a shared-scale facet_wrap). This
                # avoids the per-facet duplication that annotation_logticks
                # would produce. Subticks share the same length as decade
                # ticks (theme controls one tick length).
                .build_log_tick_ys <- function(lo, hi) {
                    out <- numeric(0)
                    k_min <- floor(log10(lo))
                    k_max <- ceiling(log10(hi))
                    for (k in k_min:k_max) {
                        for (m in 2:9) {
                            v <- m * 10^k
                            if (v > lo && v < hi) out <- c(out, v)
                        }
                    }
                    out
                }
                .log_subtick_ys <- .build_log_tick_ys(miny, axis_top)
                all_breaks_log  <- sort(unique(c(breaks_vec, .log_subtick_ys)))
                .is_major       <- all_breaks_log %in% breaks_vec
                all_labels_log  <- ifelse(.is_major,
                                          .fmt_log_break(all_breaks_log),
                                          "")

                p_data <- p_data +
                    ggplot2::scale_y_continuous(
                        trans  = "log10",
                        limits = c(miny, axis_top),
                        breaks = all_breaks_log,
                        labels = all_labels_log,
                        expand = ggplot2::expansion(mult = c(y_expand_low_mult, 0)),
                        oob    = scales::oob_keep
                    ) +
                    ggplot2::coord_cartesian(clip = "off")

                # (Subticks rendered via the scale's `breaks` argument
                # above, so the THEME's axis.ticks.y mechanism draws them
                # -- automatically respects facet behavior: only the
                # leftmost facet in a shared-scale facet_wrap shows
                # them. Trade-off: subticks share the decade ticks'
                # length since theme controls one length.)
            } else if (isTRUE(tight_top_axis)) {
                breaks_vec <- pretty(c(0, force_y_top * 1.001), n = 5)
                top_break  <- max(breaks_vec)
                axis_top   <- top_break * (1 + tight_top_axis_overhang)
                p_data <- p_data +
                    ggplot2::scale_y_continuous(
                        breaks = breaks_vec,
                        limits = c(0, axis_top),
                        expand = ggplot2::expansion(mult = c(y_expand_low_mult, 0))
                    )
            } else {
                breaks_vec <- pretty(c(0, force_y_top), n = 5)
                p_data <- p_data +
                    ggplot2::scale_y_continuous(
                        limits = c(0, NA),
                        expand = ggplot2::expansion(mult = c(y_expand_low_mult, y_expand_high_mult))
                    )
            }

            # ---- Widest y-axis label, for matching the strip-plot chrome --
            # The stats strip plot (built below) renders y-axis chrome with
            # color = NA so it's invisible but still occupies space. We need
            # the strip's invisible labels to be at least as WIDE as the
            # widest visible label on the data plot, otherwise patchwork
            # stacks the two plots with mismatched left-margin widths and
            # the brackets sit a few pixels to the LEFT of the boxplot
            # centres. By forcing the strip to render an invisible copy of
            # the data plot's widest label, both plots reserve identical
            # left-margin space and align panel-by-panel.
            .widest_y_label <- {
                labs_chr <- as.character(breaks_vec)
                # Prefer the longest character string; if breaks_vec is
                # empty (defensive), fall back to a 5-char placeholder.
                if (length(labs_chr) > 0) labs_chr[which.max(nchar(labs_chr))]
                else                       "0.001"
            }

            if (!is.null(fill_palette)) {
                p_data <- p_data + ggplot2::scale_fill_manual(
                    values = fill_palette, name = peptide_col, drop = FALSE)
            }

            # Facet with strip on bottom in "group" mode (CMV labels below
            # each panel); otherwise default top strip carries the group label.
            strip_position_data <- if (identical(x_axis_mode, "group")) "bottom" else "top"

            # ALIGNMENT FIX: force_panelsizes REMOVED from both p_data and
            # p_strip in the shared-y / build_cluster_facet path.
            #
            # WHY: when two MULTI-FACETED plots both use
            # ggh4x::force_panelsizes(cols = grid::unit(cm)) and get
            # stacked via patchwork, the resulting gtable composition
            # can apply a horizontal scale to one of the plots to make
            # totals reconcile -- the "strip plot looks compressed in x"
            # symptom. The same `force_panelsizes(cols)` works fine when
            # both stacked plots are SINGLE-PANEL (see
            # plotAbundancesDiffStrip.R) but fails on multi-faceted
            # stacks. Without force_panelsizes, patchwork's natural
            # column-width alignment kicks in: both plots have the same
            # facets, the same x scale, and the same chrome, so panel
            # columns get identical widths automatically.
            #
            # TRADE-OFF: `panel_width_cm` no longer controls panel width
            # in this mode. Panel widths now depend on the figure width
            # and chrome widths. The user controls overall size via
            # fig.width in the chunk and the n_cols / cluster_block_gap_pt
            # parameters.
            .n_facet_cols <- length(group_levels)
            p_data <- p_data +
                ggplot2::scale_x_discrete(limits = peptide_levels, drop = FALSE) +
                ggh4x::facet_wrap2(
                    facets        = ggplot2::vars(!!rlang::sym(group_col)),
                    nrow          = 1,
                    strip.position = strip_position_data
                ) +
                ggplot2::labs(x = NULL, y = y_axis_root) +
                ggplot2::theme_bw() +
                ggplot2::theme(
                    panel.grid       = ggplot2::element_blank(),
                    panel.spacing    = grid::unit(0, "lines"),
                    # aspect.ratio removed -- force_panelsizes(rows + cols)
                    # now locks panel shape absolutely.
                    plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5,
                                                             size = title_size,
                                                             margin = ggplot2::margin(b = 2)),
                    axis.title.y     = ggplot2::element_text(size = row_label_size, face = "bold"),
                    axis.text        = ggplot2::element_text(color = "black", size = axis_text_size),
                    # Push y-axis tick labels rightward (away from the
                    # axis/ticks) so they don't crowd the log subticks.
                    axis.text.y      = ggplot2::element_text(
                                            color  = "black",
                                            size   = axis_text_size,
                                            margin = ggplot2::margin(r = 4)
                                        ),
                    axis.text.x      = if (identical(x_axis_mode, "group")) {
                                            ggplot2::element_blank()
                                       } else {
                                            ggplot2::element_text(angle = 45, hjust = 1, vjust = 1,
                                                                  size = axis_text_size)
                                       },
                    # In group mode the x-axis line AND the centre tick are
                    # drawn manually via geom_segment below so they share
                    # the same y position. Blank the theme-controlled ticks
                    # to avoid a stray "second tick" at the panel bottom.
                    axis.ticks.x     = if (identical(x_axis_mode, "group"))
                                            ggplot2::element_blank()
                                       else ggplot2::element_line(color = "black"),
                    axis.ticks       = ggplot2::element_line(color = "black"),
                    # Strip: CMV-status label per facet. When at bottom (group
                    # mode), rotate 45 deg for compactness.
                    strip.background = ggplot2::element_rect(fill = NA, color = NA),
                    strip.text.x     = ggplot2::element_text(
                        size  = axis_text_size,
                        face  = "plain",
                        angle = if (identical(x_axis_mode, "group")) x_axis_label_angle else 0,
                        hjust = x_axis_label_hjust,
                        vjust = x_axis_label_vjust,
                        margin = if (identical(x_axis_mode, "group"))
                                     ggplot2::margin(t = x_axis_label_margin_top, b = 2)
                                 else ggplot2::margin(b = 2, t = 2)
                    ),
                    plot.margin = ggplot2::margin(2, 2, 2, 2)
                )

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

            # ---- Manual x-axis line inside each facet ----------------------
            # In group mode, strip.position = "bottom" places the strip text
            # BELOW the panel but ABOVE the theme-drawn axis.line.x. That
            # leaves a visible gap between the panel bottom and the axis
            # line, which then doesn't meet the y-axis at the bottom-left
            # corner. Fix: hide the theme's axis.line.x and draw a manual
            # geom_segment INSIDE each facet, at the bottom of the data
            # range.
            #
            # `x = -Inf, xend = Inf` is critical here. The discrete x scale's
            # default expansion adds 0.6 padding on each side of the data
            # range (positions 1..3 -> panel covers 0.4 to 3.6). Using
            # finite x values like 0.5 / 3.5 leaves a 0.1 gap to each panel
            # boundary; -Inf / Inf extend to the actual panel edges so
            # adjacent panels' segments butt up at the facet boundary
            # (panel.spacing = 0) for a genuinely continuous axis line.
            if (identical(x_axis_mode, "group")) {
                # y position of the axis line: NUDGED slightly below the
                # minimum data position so points at miny (log) / 0 (linear)
                # sit clearly above the line instead of crossing it.
                #
                # CRITICAL: the offset must be a FRACTION of the panel's
                # lower-expansion margin (NOT a fixed log-unit / linear-
                # unit offset). Otherwise short-y-range clusters (e.g. a
                # low-abundance combined_merg cluster with y data in
                # [0.1, 0.5]) push the manual axis line + tick PAST the
                # panel bottom -- they then render outside the panel via
                # coord_cartesian(clip = "off"), leaving a visible gap
                # between the y-axis and the (orphaned) x-axis line, and
                # making the tick look ridiculously long. Setting the
                # offsets as fractions of the lower expansion margin
                # makes the visual offset / tick length identical in
                # PIXELS regardless of the y data scale.
                #
                # Axis line at 75% of the lower expansion margin
                # (closer to panel bottom than halfway -- gives points
                # sitting at miny / 0 more vertical clearance above the
                # axis line so they don't visually clip it).
                # Tick extends to 95% of the margin (almost at panel
                # bottom). Both stay safely INSIDE the panel.
                .log_range <- if (isTRUE(log)) log10(axis_top) - log10(miny) else NA_real_
                axis_y_value <- if (isTRUE(log)) {
                                    miny * 10 ^ (-0.75 * y_expand_low_mult * .log_range)
                                } else {
                                    -force_y_top * y_expand_low_mult * 0.75
                                }
                # The centre tick extends a bit further below the line.
                tick_y_end   <- if (isTRUE(log)) {
                                    miny * 10 ^ (-0.95 * y_expand_low_mult * .log_range)
                                } else {
                                    -force_y_top * y_expand_low_mult * 0.95
                                }

                n_groups_local2 <- length(group_levels)
                axis_line_df <- data.frame(
                    x    = rep(-Inf, n_groups_local2),
                    xend = rep( Inf, n_groups_local2),
                    y    = rep(axis_y_value, n_groups_local2),
                    yend = rep(axis_y_value, n_groups_local2),
                    stringsAsFactors = FALSE
                )
                axis_line_df[[group_col]] <- factor(group_levels, levels = group_levels)

                # Centre tick per facet: one tick at the middle x position
                # of each CMV condition, drawn as a manual segment from the
                # axis line downward. Using a manual segment (rather than
                # theme-controlled axis.ticks.x with a scale_x_discrete
                # `breaks` override) guarantees the tick connects to our
                # manual axis line at exactly the same y, with no offset
                # between them.
                center_x  <- (length(peptide_levels) + 1L) / 2
                tick_df   <- data.frame(
                    x    = rep(center_x, n_groups_local2),
                    xend = rep(center_x, n_groups_local2),
                    y    = rep(axis_y_value, n_groups_local2),
                    yend = rep(tick_y_end,   n_groups_local2),
                    stringsAsFactors = FALSE
                )
                tick_df[[group_col]] <- factor(group_levels, levels = group_levels)

                # Manual Y-axis line on the LEFTMOST facet, terminating
                # exactly at the manual x-axis line (y = axis_y_value).
                # The theme-drawn axis.line.y runs the full panel height
                # (top to bottom), so when our manual x-axis sits ABOVE
                # the panel bottom (50% into the lower expansion margin),
                # the theme y-axis line overshoots downward past the
                # x-axis line -> visible overhang at the bottom-left
                # corner. Drawing y-axis manually here lets us cap it at
                # exactly axis_y_value so the corner closes cleanly.
                #
                # `x = -Inf` puts the segment at the panel's left edge
                # (same as where theme's axis.line.y draws). Only the
                # FIRST group_col level is tagged so this segment renders
                # ONLY on the leftmost facet (matches shared_y_per_row's
                # convention of leftmost-only y-axis chrome).
                y_axis_line_df <- data.frame(
                    x    = -Inf,
                    xend = -Inf,
                    y    = axis_y_value,
                    yend = axis_top,
                    stringsAsFactors = FALSE
                )
                y_axis_line_df[[group_col]] <- factor(group_levels[1],
                                                       levels = group_levels)

                p_data <- p_data +
                    ggplot2::geom_segment(
                        data    = axis_line_df,
                        mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                        inherit.aes = FALSE,
                        color = "black", linewidth = 0.4
                    ) +
                    ggplot2::geom_segment(
                        data    = tick_df,
                        mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                        inherit.aes = FALSE,
                        color = "black", linewidth = 0.4
                    ) +
                    ggplot2::geom_segment(
                        data    = y_axis_line_df,
                        mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                        inherit.aes = FALSE,
                        color = "black", linewidth = 0.4
                    ) +
                    # Suppress BOTH theme-drawn axis lines now -- manual
                    # geom_segments take their place and meet cleanly at
                    # the bottom-left corner (no overhang).
                    ggplot2::theme(
                        axis.line.x = ggplot2::element_blank(),
                        axis.line.y = ggplot2::element_blank()
                    )

                # ---- CMV labels: auto-anchored to centre tick ------------
                # The strip-text route (theme strip.text.x) anchors the
                # label horizontally WITHIN the panel-wide strip area, so
                # hjust controls left/center/right of the panel -- you
                # cannot make the rightmost letter sit at the centre tick
                # (panel-middle) via hjust alone.
                #
                # Instead we draw the labels as a panel-scoped geom_text
                # at the centre-tick x position. With:
                #   hjust = 1, vjust = 1, angle = 45
                # the rightmost letter sits AT the anchor (which we place
                # at the tick endpoint), and the rest of the label trails
                # diagonally down-and-to-the-left -- the conventional
                # tilted-axis-label look. Strip-text rendering is then
                # suppressed (element_blank) to avoid double labels.
                #
                # Label y is placed just below the tick endpoint. The
                # user's x_axis_label_hjust / vjust / angle parameters
                # still apply (default to 1/1/45). x_axis_label_margin_top
                # is INTENTIONALLY UNUSED in this code path -- vertical
                # position is computed from axis_y_value / tick_y_end.
                .label_y_offset_mult <- if (isTRUE(log)) 10^(-0.04) else 1.05
                .label_y <- if (isTRUE(log)) tick_y_end * .label_y_offset_mult
                            else            tick_y_end - force_y_top * 0.008
                label_df <- data.frame(
                    x = rep(center_x, n_groups_local2),
                    y = rep(.label_y,  n_groups_local2),
                    label = as.character(group_levels),
                    stringsAsFactors = FALSE
                )
                label_df[[group_col]] <- factor(group_levels, levels = group_levels)

                p_data <- p_data +
                    ggplot2::geom_text(
                        data        = label_df,
                        mapping     = ggplot2::aes(x = x, y = y, label = label),
                        inherit.aes = FALSE,
                        angle       = x_axis_label_angle,
                        hjust       = x_axis_label_hjust,
                        vjust       = x_axis_label_vjust,
                        # convert pt -> ggplot mm-pt for geom_text:
                        size        = axis_text_size / 2.845
                    ) +
                    # KEEP strip.text.x as an INVISIBLE element_text so the
                    # strip row in the gtable still reserves vertical
                    # space below the panel. coord_cartesian(clip = "off")
                    # lets the geom_text labels (drawn in panel coords)
                    # extend down past the panel boundary into that
                    # reserved strip-row space. If we set strip.text.x =
                    # element_blank(), the strip row would COLLAPSE and
                    # there would be no plot-region space for the labels
                    # to render -- they'd be clipped by the plot-region
                    # boundary (since clip = "off" only disables the
                    # panel-level clip, not the plot-region clip).
                    #
                    # We use the longest plausible label ("Persistent IgM"
                    # is 14 chars at 17pt rotated 45deg => ~80-100pt of
                    # vertical strip height needed) as the placeholder so
                    # space is reserved correctly.
                    ggplot2::theme(
                        strip.background = ggplot2::element_rect(fill = NA, color = NA),
                        strip.text.x     = ggplot2::element_text(
                            size   = axis_text_size,
                            angle  = x_axis_label_angle,
                            hjust  = x_axis_label_hjust,
                            vjust  = x_axis_label_vjust,
                            color  = NA,  # invisible -- geom_text renders the actual labels
                            margin = ggplot2::margin(t = x_axis_label_margin_top, b = 2)
                        )
                    )
            }

            # ---- Stats strip ggplot with facet_wrap2 -----------------------
            # Built only if we have brackets; otherwise we emit an empty
            # ggplot of the right facet shape so heights line up.
            # Match the data plot's strip.position so the two stacked
            # gtables have IDENTICAL row layouts (strip row in same slot,
            # panel row in same slot). Patchwork aligns columns within
            # matching rows; mismatched row structures can produce small
            # horizontal offsets.
            n_x_strip <- length(peptide_levels)
            p_strip   <- ggplot2::ggplot() +
                ggh4x::facet_wrap2(
                    facets         = ggplot2::vars(!!rlang::sym(group_col)),
                    nrow           = 1,
                    strip.position = strip_position_data
                )
            if (!is.null(seg_df_all) && nrow(seg_df_all) > 0) {
                p_strip <- p_strip + ggplot2::geom_segment(
                    data    = seg_df_all,
                    mapping = ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                    color   = bracket_color, linewidth = bracket_size,
                    lineend = "round", linejoin = "round"
                )
                if (isTRUE(use_richtext)) {
                    p_strip <- p_strip + ggtext::geom_richtext(
                        data    = text_df_all,
                        mapping = ggplot2::aes(x = x, y = y, label = label, vjust = vj),
                        size    = label_size, fontface = label_fontface,
                        label.size = NA, fill = NA,
                        label.padding = grid::unit(c(0, 0, 0, 0), "lines")
                    )
                } else {
                    p_strip <- p_strip + ggplot2::geom_text(
                        data    = text_df_all,
                        mapping = ggplot2::aes(x = x, y = y, label = label, vjust = vj),
                        size    = label_size, fontface = label_fontface
                    )
                }
            }
            # Force facet panels to share the SAME group_col domain even when
            # a panel has no brackets (so panels stay aligned with data above).
            # We add a tiny invisible dummy point per group to ensure all
            # facets are realized. The data.frame is constructed with one
            # row per group_level so the factor column assignment matches.
            n_groups_local <- length(group_levels)
            dummy_strip <- data.frame(
                x = rep(1, n_groups_local),
                y = rep(0, n_groups_local)
            )
            dummy_strip[[group_col]] <- factor(group_levels, levels = group_levels)
            p_strip <- p_strip + ggplot2::geom_blank(
                data    = dummy_strip,
                mapping = ggplot2::aes(x = x, y = y),
                inherit.aes = FALSE
            )

            # TIGHT y-limits so brackets fill the strip's allocated
            # panel height, eliminating empty whitespace between the
            # bottommost bracket and the data plot's panel below.
            #
            # Bracket geometry (in strip y units):
            #   - bracket bar at y = bracket_y[k] = (n - k + 0.5)
            #   - bracket tip extends DOWN by bracket_tip_npc
            #   - bracket label extends UP by label_nudge_npc (+ asterisk
            #     offset for *-labels)
            # So the actual y range used is roughly:
            #   [0.5 - bracket_tip_npc, n_max - 0.5 + label_nudge_npc + ~0.35]
            # We pad each end by a tiny amount so tips/labels aren't clipped.
            .strip_y_lo <- (0.5 - bracket_tip_npc) - 0.02
            .strip_y_hi <- (y_axis_max_strip - 0.5) +
                            label_nudge_npc + asterisk_y_offset + 0.35
            # Defensive: if no brackets in this cluster, the formulas
            # collapse to a tiny range; widen it slightly so the panel
            # still renders.
            if (max_brackets_cl == 0L) {
                .strip_y_lo <- 0
                .strip_y_hi <- 1
            }

            p_strip <- p_strip +
                # CRITICAL: use IDENTICAL scale type to the data plot
                # (scale_x_discrete with the SAME limits). Trying to fake
                # a matching panel x range via scale_x_continuous +
                # bespoke expansion is fragile -- discrete and continuous
                # scales compute panel ranges and break-positions
                # differently in subtle ways (e.g. how expand interacts
                # with the discrete data range, how mapped_discrete
                # numerics get positioned). Using the SAME discrete scale
                # on both plots guarantees position k maps to the same
                # panel fraction in both plots, so brackets land exactly
                # over boxplot centres.
                #
                # We pass numeric x values (integers and midpoints like
                # 1.5) to a discrete scale: ggplot2 >= 3.3 handles this
                # via `mapped_discrete` -- numerics are interpreted as
                # positions on the discrete axis (1 -> first level,
                # 1.5 -> halfway between first and second, etc.).
                ggplot2::scale_x_discrete(
                    limits = peptide_levels, drop = FALSE
                ) +
                ggplot2::scale_y_continuous(
                    limits = c(.strip_y_lo, .strip_y_hi),
                    # Force the strip plot to render an invisible y-axis
                    # label equal in width to the data plot's widest visible
                    # label. axis.text.y has color = NA below, so the label
                    # itself is invisible; only its bounding box (and the
                    # left-margin space it reserves) is preserved. This is
                    # what gets patchwork's column alignment to match the
                    # data plot's leftmost panel position exactly so the
                    # brackets sit centered over the boxplots below.
                    breaks = c((.strip_y_lo + .strip_y_hi) / 2),
                    labels = c(.widest_y_label),
                    expand = c(0, 0), oob = scales::oob_keep
                ) +
                ggplot2::coord_cartesian(clip = "off") +
                # MATCH the data plot's labs(x = NULL, y = ...) so the strip
                # plot also collapses the axis.title.x row and reserves a
                # y-title row of the same nominal size.
                ggplot2::labs(x = NULL, y = y_axis_root) +
                # force_panelsizes REMOVED -- see matching comment on the
                # data plot above. With both plots free of force_panelsizes
                # in this multi-faceted stack, patchwork's natural column
                # alignment makes their panel widths match exactly.
                # (panel_width_cm is therefore unused in this mode.)
                # IMPORTANT: don't use theme_void() here. The data plot has
                # y-axis chrome (axis title + tick labels + tick marks + axis
                # line) on its leftmost facet, which takes up horizontal
                # space to the LEFT of its leftmost panel. If the strip
                # plot lacks that chrome, patchwork stacks the two plots
                # with the strip's leftmost panel hugging the left of the
                # plot region, while the data plot's leftmost panel is
                # pushed right by its chrome width -> the strip's brackets
                # appear shifted LEFT relative to the boxplots below.
                #
                # Fix: reserve the same y-axis chrome (invisibly) so the
                # strip's leftmost panel starts at the same x as the data
                # plot's leftmost panel. We achieve this with explicit
                # element_text(color = NA) for axis.title.y and axis.text.y
                # and element_line(color = NA) for axis.ticks.y / axis.line.y.
                # (labs already set above to match data plot's labs)
                ggplot2::theme_bw() +
                ggplot2::theme(
                    # --- y-axis chrome: VISIBLE on data plot -> invisible
                    #     element_text on strip so chrome WIDTH still matches
                    #     (panel-column alignment between the two plots). -----
                    axis.title.y     = ggplot2::element_text(
                                            size  = row_label_size,
                                            face  = "bold",
                                            color = NA
                                        ),
                    axis.text.y      = ggplot2::element_text(
                                            color  = NA,
                                            size   = axis_text_size,
                                            margin = ggplot2::margin(r = 4)
                                        ),
                    axis.ticks.y     = ggplot2::element_line(color = NA),
                    axis.line.y      = ggplot2::element_line(color = NA, linewidth = 0.4),
                    # --- x-axis chrome and strip row: ALL element_blank() so
                    #     the corresponding gtable rows COLLAPSE entirely.
                    #     This squeezes out the vertical whitespace between
                    #     the brackets and the data plot's panel below.
                    #     We no longer need to reserve these rows for gtable-
                    #     structure matching now that force_panelsizes is off
                    #     and patchwork's natural column alignment handles x.
                    axis.title.x     = ggplot2::element_blank(),
                    axis.text.x      = ggplot2::element_blank(),
                    axis.ticks.x     = ggplot2::element_blank(),
                    axis.line.x      = ggplot2::element_blank(),
                    strip.background = ggplot2::element_blank(),
                    strip.text.x     = ggplot2::element_blank(),
                    # --- panel chrome (no border, no grid, no bg) -----------
                    panel.border     = ggplot2::element_blank(),
                    panel.background = ggplot2::element_blank(),
                    panel.grid       = ggplot2::element_blank(),
                    plot.background  = ggplot2::element_blank(),
                    # aspect.ratio dropped: force_panelsizes(rows + cols)
                    # above locks panel shape absolutely so it can't get
                    # compressed by patchwork's vertical-space budget.
                    panel.spacing    = grid::unit(0, "lines"),
                    plot.margin      = ggplot2::margin(0, 2, 0, 2)
                )

            # ---- Compose strip / data for this cluster ---------------------
            cluster_comp <- p_strip / p_data +
                patchwork::plot_layout(heights = c(strip_h, data_h))

            cluster_comp
        }  # end build_cluster_facet()

        # Build per-cluster blocks (composite + header).
        # `is_leftmost_in_row` controls whether THIS cluster block shows
        # the y-axis title ("Proportion [%]" by default). Only the
        # leftmost cluster in each row gets it; others have it blanked
        # to avoid repeating the same label across the row.
        cluster_blocks <- list()
        for (cl_i in seq_along(clusters)) {
            cl <- clusters[cl_i]
            cluster_comp <- build_cluster_facet(cl)

            # Determine row position. For "wrap" layout, blocks fill
            # row-by-row at n_cols per row; for "single_row" everything
            # is in one row so only the first cluster is leftmost.
            .is_leftmost_in_row <- if (identical(layout, "wrap")) {
                ((cl_i - 1L) %% n_cols) == 0L
            } else {  # "single_row"
                cl_i == 1L
            }
            if (!.is_leftmost_in_row) {
                cluster_comp <- cluster_comp + ggplot2::theme(
                    axis.title.y = ggplot2::element_blank()
                )
            }

            header_p <- ggplot2::ggplot() +
                ggplot2::theme_void() +
                ggplot2::labs(title = cl) +
                ggplot2::theme(
                    plot.title = ggplot2::element_text(
                        hjust  = 0.5,
                        face   = cluster_header_face,
                        size   = cluster_header_size,
                        margin = ggplot2::margin(t = 2, b = 2)
                    )
                )

            block <- header_p / cluster_comp +
                patchwork::plot_layout(
                    heights = grid::unit.c(grid::unit(1.5, "lines"),
                                           grid::unit(1, "null"))
                ) &
                ggplot2::theme(
                    plot.margin = ggplot2::margin(
                        t = 0, r = cluster_block_gap_pt,
                        b = 0, l = if (cl_i == 1L) 0 else cluster_block_gap_pt,
                        unit = "pt"
                    )
                )

            cluster_blocks[[cl]] <- block
        }

        # Wrap into the final figure:
        #   layout = "single_row" -> nrow = 1 (all blocks in one row).
        #     `n_cols` is IGNORED on this path -- prevents the empty-slot
        #     artifact where wrap_plots(ncol = N > n_clusters) reserves
        #     N column widths and the legend overlaps the empty slot.
        #   layout = "wrap"       -> ncol = n_cols (wraps to multiple
        #     rows). Use for many-cluster figures.
        out <- if (identical(layout, "single_row")) {
            patchwork::wrap_plots(cluster_blocks, nrow = 1)
        } else {  # "wrap"
            patchwork::wrap_plots(cluster_blocks, ncol = n_cols)
        }
        out <- out &
            .build_legend_theme(legend_position,
                                legend_key_size, legend_text_size,
                                legend_title_size, legend_spacing_y_pt,
                                axis_text_size, row_label_size) &
            .build_fill_guide(peptide_col, legend_key_size)
        if (legend_position != "none") {
            out <- out + patchwork::plot_layout(guides = "collect")
        }

    } else {
        # ---- Original grid layout (rows = clusters, cols = groups) -----------
        for (cl in clusters) for (gl in group_levels) {
            idx <- idx + 1L
            col_pos <- ((idx - 1L) %% n_cols_grid) + 1L
            composites[[idx]] <- make_composite(
                cl, gl, col_pos = col_pos,
                include_cluster_in_y_label = TRUE,
                show_col_title_for_strip   = (((idx - 1L) %/% n_cols_grid) + 1L) == 1L
            )
        }

        out <- patchwork::wrap_plots(composites, ncol = n_cols_grid) &
            .build_legend_theme(legend_position,
                                legend_key_size, legend_text_size,
                                legend_title_size, legend_spacing_y_pt,
                                axis_text_size, row_label_size) &
            .build_fill_guide(peptide_col, legend_key_size)
        if (legend_position != "none") {
            out <- out + patchwork::plot_layout(guides = "collect")
        }
    }

    # ---- Methods-explainer title (toggleable) ---------------------------
    # patchwork-level annotation describing the stats; sits above the whole
    # composite. Plays nicely with f2() — the title is part of the saved PDF.
    if (!is.null(stats_title) && !identical(stats_title, FALSE)) {
        # Derive the correction label from the actively-used label_col so
        # the title accurately reflects what the stars represent. The method
        # name (Holm / BH / etc.) comes from the function's p_adjust_method.
        .method_label <- switch(p_adjust_method,
            holm       = "Holm-Bonferroni correction",
            BH         = "Benjamini-Hochberg (BH) FDR correction",
            bonferroni = "Bonferroni correction",
            hochberg   = "Hochberg correction",
            hommel     = "Hommel correction",
            BY         = "Benjamini-Yekutieli (BY) FDR correction",
            paste0("'", p_adjust_method, "' correction")
        )
        .correction_label <- if (identical(label_col, "p.adj_global.signif")) {
            paste0("global ", .method_label, " (across the full figure)")
        } else if (identical(label_col, "p.adj.signif") ||
                   identical(label_col, "p.adj_per_cluster.signif")) {
            paste0("per-cluster ", .method_label)
        } else {
            paste0("uses column '", label_col, "' as the significance source")
        }

        title_text <- if (isTRUE(stats_title)) {
            paste0(
                "Paired Wilcoxon signed-rank tests vs '", reference_level,
                "' per donor; ", .correction_label,
                ".  Stars: * p<0.05  ** p<0.01  *** p<0.001  **** p<0.0001",
                if (isTRUE(hide_ns)) "  (ns hidden)" else ""
            )
        } else {
            as.character(stats_title)
        }

        out <- out + patchwork::plot_annotation(
            title = title_text,
            theme = ggplot2::theme(
                plot.title = ggplot2::element_text(
                    size   = stats_title_size,
                    face   = stats_title_face,
                    hjust  = 0,
                    margin = ggplot2::margin(b = 6)
                )
            )
        )
    }

    # ---- Stash source data + full stats on the patchwork ----------------
    attr(out, "source_data")  <- df
    attr(out, "source_stats") <- full_stats

    if (isTRUE(verbose)) {
        axis_width_cm   <- 1.4
        legend_width_cm <- if (legend_position == "right") 3.0 else 0
        gap_cm          <- panel_gap_pt * 0.03528
        col_width_cm    <- panel_width_cm + 2 * gap_cm
        total_cm        <- n_cols_grid * col_width_cm + axis_width_cm + legend_width_cm
        total_in        <- total_cm / 2.54
        message(sprintf(
            "plotAimPaired: estimated min fig.width = %.1f in (%.1f cm) [n_cols=%d, panel_width_cm=%g, legend='%s']. Increase if PDFs look clipped.",
            total_in, total_cm, n_cols_grid, panel_width_cm, legend_position
        ))
    }

    out
}

# Attach CATALYST namespace for .check_k / cluster_ids resolution.
