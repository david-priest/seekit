library(SingleCellExperiment)
library(ggplot2)
library(ggbeeswarm)
library(scales)
library(dplyr)
library(rstatix)
library(rlang)
library(ggh4x)
library(ggnewscale)

# Plot pseudobulk expression split by a cell-level activation clustering (act_k),
# faceted by a cell-level phenotype clustering (pheno_k).
#
# Each point = median expression for one sample within a (pheno_cluster, act_cluster) bin.
# X-axis    = act_k levels
# Fill      = act_k levels
# Facet rows = pheno_k clusters
# Facet cols = antigens (features)
#
# Both pheno_k and act_k are resolved from colData(x) directly, or from
# .wl_cluster_codes(x) if not found in colData.

plotPbActExprsDiff <- function(
    x,
    pheno_k        = "merging_pheno",  # cell-level colData col for facet rows
    act_k          = "merging_act",    # cell-level colData col for x-axis groups
    features       = "state",
    assay          = "exprs",
    fun            = c("median", "mean", "sum"),
    conditions     = NULL,             # filter act_k levels to include
    excluded_pheno = NULL,             # exclude pheno cluster levels
    excluded_act   = NULL,             # exclude act cluster levels
    point_size     = 1.5,
    textsize       = 14,
    panel_spacing  = 2,
    show_stats     = FALSE,
    hide_ns        = TRUE,
    stat_size      = 4,
    nudge          = 0,
    step_increase  = 0.3,
    group1keep     = NULL,
    external_stats = NULL,
    facet_ratio    = 1.5,
    scales         = "free_y",
    geom           = c("boxes", "bar"),
    swap           = FALSE,            # swap facet rows/cols
    hide_x_labels  = FALSE,
    # When TRUE, draws each pheno_cluster name as a 45-deg rotated
    # geom_text label below the panel, with the rightmost letter
    # anchored AT the centre x-axis tick (single tick per panel at the
    # middle act_cluster position) -- same auto-anchoring trick as
    # plotAimPaired's group mode. The per-panel act-cluster tick labels
    # are hidden, the strip text at the column bottom is suppressed,
    # `axes` is forced to "x" so the x-axis line is drawn on every
    # antigen row, panel.spacing.x is set to 0 so adjacent panels sit
    # flush, and the boxplots within each panel are widened to fill the
    # panel (controlled via `box_width_grouped` below). Only meaningful
    # when swap = FALSE (pheno_cluster on cols).
    cluster_label_below = FALSE,
    # Boxplot width used when cluster_label_below = TRUE (boxes within a
    # panel sit tightly together for the "panel = one cluster's three
    # bars" look). Default 0.85 gives slight breathing room between the
    # 3 boxes in each panel; bump to 0.95 to make them nearly abut.
    box_width_grouped   = 0.85,
    # Horizontal spacing BETWEEN adjacent pheno_cluster panels (the
    # "sets of 3 boxplots") when cluster_label_below = TRUE. pt units.
    # 0 makes adjacent panels abut and the x-axis line look continuous;
    # bump to ~6-10 pt for visible separation between clusters.
    # Note: any value > 0 introduces a small gap in the per-panel x-axis
    # line between panels (since each panel draws its own line).
    cluster_panel_spacing_pt = 8,
    # Vertical spacing BETWEEN antigen panel rows when cluster_label_below
    # = TRUE. pt units. NULL (default) preserves the generic `panel_spacing`
    # (lines) value already set in the base theme block. Override with an
    # explicit pt value when point glyphs at the top of one panel are
    # spilling into the panel above -- the tight_top_axis termination
    # leaves only ~2% headroom past data_max, so a finite point radius can
    # render above the axis line and visually reach the next row.
    cluster_panel_spacing_y_pt = NULL,
    # Auto-anchor controls for the rotated cluster labels below each
    # panel. Defaults reproduce plotAimPaired's "right letter at the
    # centre tick" look.
    cluster_label_angle = 45,
    cluster_label_hjust = 1,
    # vjust slightly > 1 shifts the rotated label DOWN from the axis by
    # (vjust - 1) * bbox_height. 1.2 gives ~3-5 pixels of breathing room
    # between the axis line and the top edge of the rotated text at
    # ~17-18pt sizes. Bump to e.g. 1.4 for more separation.
    cluster_label_vjust = 1.2,
    cluster_label_size  = NULL,  # NULL -> textsize
    # ---- cluster separators (dashed verticals between pheno_cluster groups) ----
    # Continuous dashed vertical line in EACH gutter between adjacent
    # pheno_cluster panel columns, spanning the FULL panel-area height
    # (all antigen rows + the panel.spacing.y rows between them). Only
    # meaningful when there are >= 2 panel columns; defaults to off so
    # existing callers don't change behaviour.
    #
    # Implementation: a geom_vline inside the ggplot can only span one
    # panel at a time, so the inter-row panel.spacing.y gutters would
    # break the visual line. We instead overlay a grid::segmentsGrob at
    # the gtable level (in the gap layout-column between adjacent panel
    # columns, spanning multiple gtable rows) -- but to keep the function
    # returning a ggplot the user can still `+ scale_fill_manual(...)` to,
    # we DEFER the gtable surgery until print time via an S3 print method
    # (print.plotPbActExprsDiff, defined at the bottom of this file).
    # When cluster_separator = TRUE, the returned plot has class
    # `c("plotPbActExprsDiff", "gg", "ggplot")` and the separator config
    # stashed as attribute `.cluster_separator_cfg`.
    cluster_separator          = FALSE,
    cluster_separator_linetype = "dashed",
    cluster_separator_color    = "grey50",
    # grid lwd units (~ points). 1.0 is roughly half a typical ggplot
    # axis-line width -- visible but unobtrusive. Bump if dashes don't
    # render on low-DPI previews.
    cluster_separator_lwd      = 1.0,
    label_parse    = FALSE,
    axes           = "all",
    point_color_by = NULL,             # sample-level colData col to colour points by
    dodge          = FALSE,            # if TRUE and point_color_by is set, dodge boxes and
                                       # points by point_color_by within each act_cluster group
    pheno_order    = NULL,             # character vector: desired order of pheno_cluster facets
    act_order      = NULL,             # character vector: desired order of act_cluster x-axis groups
    # ---- connect-donors / paired-slope overlay -------------------------------
    # When non-NULL, draws grey lines connecting same-donor points across the
    # act_cluster x positions within each (pheno_cluster, antigen) facet cell.
    # The named column must be sample-level (constant within sample_id, e.g.
    # "patient_id" or "sample_id").
    #
    # Implementation follows the plotAbundancesA1 pattern (see ledger C7):
    #   1. Precompute a per-row x position .x_pos = factor(act_cluster) + jitter
    #      and bake it into the data frame.
    #   2. Replace geom_quasirandom with geom_point(aes(x = .x_pos)) so the
    #      points sit at the same x's the line endpoints read.
    #   3. Add geom_line(aes(x = .x_pos, group = interaction(pheno, antigen,
    #      donor))) BEFORE the points so points render on top of lines.
    # NEVER use position_dodge on the geom_line layer -- position_dodge derives
    # its slot from the group aesthetic, and aes(group = donor) collapses lines
    # to vertical. Bake x into the data instead.
    #
    # Incompatible with dodge=TRUE (which uses position_dodge for points).
    # Incompatible with point_color_by != NULL in non-dodged mode (point fill
    # would still work but the layout doesn't get extra value from lines on
    # top of a multi-coloured cohort split).
    connect_donors    = NULL,
    connect_color     = "grey60",
    connect_alpha     = 0.5,
    connect_linewidth = 0.3,
    # Deterministic within-act-cluster jitter so points don't fully overlap.
    # 0 disables jitter (points sit dead-centre on each act_cluster category).
    connect_jitter_width = 0.1,
    connect_jitter_seed  = 42,
    # ---- 2026 visual-consistency parameters ----
    nature_style            = TRUE,    # blank panel border + axis lines only (matches strip plotters)
    tight_top_axis          = TRUE,    # terminate each (antigen × pheno) facet at next pretty tick above data max
    tight_top_axis_overhang = 0.02     # multiplicative overhang past topmost tick
) {
    fun  <- match.arg(fun)
    geom <- match.arg(geom)

    .wl_check_assay(x, assay)

    # ---- resolve pheno clustering ----------------------------------------
    if (pheno_k %in% names(colData(x))) {
        pheno_ids <- factor(colData(x)[[pheno_k]])
    } else {
        .wl_check_sce(x)
        pheno_ids <- .wl_cluster_ids(x, .wl_check_k(x, pheno_k))
    }

    # ---- resolve act clustering ------------------------------------------
    if (act_k %in% names(colData(x))) {
        act_ids <- factor(colData(x)[[act_k]])
    } else {
        .wl_check_sce(x)
        act_ids <- .wl_cluster_ids(x, .wl_check_k(x, act_k))
    }

    # ---- subset to requested features ------------------------------------
    x <- x[.wl_get_features(x, features), ]

    expr_mat <- as.matrix(assay(x, assay))  # features x cells

    # ---- build tidy per-cell data frame ----------------------------------
    cell_df <- data.frame(
        sample_id     = x$sample_id,
        pheno_cluster = pheno_ids,
        act_cluster   = act_ids,
        stringsAsFactors = FALSE
    )

    # Include point_color_by variable if provided (sample-level, so constant within group)
    if (!is.null(point_color_by)) {
        stopifnot(point_color_by %in% names(colData(x)))
        cell_df[[point_color_by]] <- colData(x)[[point_color_by]]
    }

    # Include connect_donors variable (sample-level; needed for per-donor line
    # groups in the (pheno_cluster × antigen) facet cells).
    if (!is.null(connect_donors)) {
        if (!(connect_donors %in% names(colData(x)))) {
            warning("plotPbActExprsDiff: connect_donors column '", connect_donors,
                    "' not found in colData; donor lines skipped.", call. = FALSE)
            connect_donors <- NULL
        } else if (connect_donors != "sample_id") {
            cell_df[[connect_donors]] <- colData(x)[[connect_donors]]
        }
    }

    agg_fun <- match.fun(fun)

    df_list <- lapply(rownames(expr_mat), function(feat) {
        tmp         <- cell_df
        tmp$value   <- expr_mat[feat, ]
        tmp$antigen <- feat
        tmp
    })
    df_full <- do.call(rbind, df_list)

    # ---- aggregate: one row per (antigen, pheno_cluster, act_cluster, sample_id)
    # connect_donors is sample-level so adding it to group_vars carries it
    # through summarise without changing the per-row identity (still 1 row
    # per sample within each panel cell).
    group_vars <- c("antigen", "pheno_cluster", "act_cluster", "sample_id",
                    point_color_by,
                    if (!is.null(connect_donors) && connect_donors != "sample_id") connect_donors)
    df <- df_full %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
        dplyr::summarise(value = agg_fun(value), n_cells = dplyr::n(), .groups = "drop") %>%
        dplyr::filter(n_cells > 0)

    DFCOND <<- df  # expose for inspection

    # ---- filter / exclude ------------------------------------------------
    if (!is.null(conditions)) {
        df <- df[df$act_cluster %in% conditions, ]
        df$act_cluster <- factor(df$act_cluster, levels = conditions, ordered = TRUE)
    }
    if (!is.null(excluded_pheno)) df <- df[!df$pheno_cluster %in% excluded_pheno, ]
    if (!is.null(excluded_act))   df <- df[!df$act_cluster   %in% excluded_act,   ]

    # ---- cluster ordering ------------------------------------------------
    if (!is.null(act_order)) {
        missing_act <- setdiff(act_order, unique(df$act_cluster))
        if (length(missing_act) > 0)
            warning("act_order contains levels not present in data: ",
                    paste(missing_act, collapse = ", "))
        df$act_cluster <- factor(df$act_cluster, levels = act_order)
    } else {
        df$act_cluster <- factor(df$act_cluster)
    }

    if (!is.null(pheno_order)) {
        missing_pheno <- setdiff(pheno_order, unique(df$pheno_cluster))
        if (length(missing_pheno) > 0)
            warning("pheno_order contains levels not present in data: ",
                    paste(missing_pheno, collapse = ", "))
        df$pheno_cluster <- factor(df$pheno_cluster, levels = pheno_order)
    } else {
        df$pheno_cluster <- factor(df$pheno_cluster)
    }

    # ---- n samples per condition (Immunity figure legend reporting) -------
    # Strategy:
    #   - If `point_color_by` is set, it's a sample-level grouping variable
    #     -> report n distinct samples per point_color_by level (the standard
    #     "patient cohort condition" report).
    #   - Otherwise fall back to total distinct sample count. `act_cluster`
    #     is a cell-level grouping (every sample contributes cells to
    #     potentially every act_cluster) so it isn't useful as a
    #     "samples-per-condition" denominator.
    if (!is.null(point_color_by) && point_color_by %in% names(df)) {
        n_per_cond_df <- df %>%
            dplyr::distinct(sample_id, !!rlang::sym(point_color_by)) %>%
            dplyr::count(!!rlang::sym(point_color_by), name = "n_samples") %>%
            as.data.frame()
        message("plotPbActExprsDiff: n samples per ", point_color_by, ":")
        print(n_per_cond_df, row.names = FALSE)
    } else {
        n_total <- length(unique(df$sample_id))
        n_per_cond_df <- data.frame(group = "all", n_samples = n_total)
        message("plotPbActExprsDiff: total n samples = ", n_total,
                " (no point_color_by set; act_cluster is cell-level, not per-sample)")
    }

    # ---- stats -----------------------------------------------------------
    if (show_stats) {
        if (!is.null(external_stats)) {
            external_stats <- external_stats %>% dplyr::filter(antigen %in% unique(df$antigen))
            dummy_stat <- df %>%
                dplyr::group_by(antigen, pheno_cluster) %>%
                rstatix::wilcox_test(value ~ act_cluster) %>%
                rstatix::add_y_position(scales = "free", step.increase = step_increase)
            external_stats <- external_stats %>%
                dplyr::left_join(
                    dummy_stat %>% dplyr::select(antigen, pheno_cluster, group1, group2, y.position),
                    by = c("antigen", "pheno_cluster", "group1", "group2")
                )
            stat.test <- external_stats
        } else {
            stat.test <- df %>%
                dplyr::group_by(antigen, pheno_cluster) %>%
                rstatix::wilcox_test(value ~ act_cluster) %>%
                rstatix::add_y_position(scales = "free", step.increase = step_increase)
        }

        if (!is.null(group1keep)) {
            stat.test <- stat.test[stat.test$group1 %in% group1keep, ]
            stat.test <- stat.test %>% rstatix::add_y_position(scales = "free", step.increase = step_increase)
        }

        pbCondStats <<- stat.test  # expose for inspection
    }

    # ---- plot ------------------------------------------------------------
    y_label <- paste(fun, ifelse(assay == "exprs", "expression", assay))

    p <- ggplot(df, aes(x = act_cluster, y = value)) +
        labs(x = act_k, y = NULL) +  # y-axis title suppressed; antigen strip on left acts as label
        theme_bw() +
        theme(
            panel.grid        = element_blank(),
            text              = element_text(size = textsize),
            strip.text        = element_text(size = textsize),
            strip.text.y.left = element_text(size = textsize, angle = 90, hjust = 0.5),
            strip.background  = element_rect(fill = NA, color = NA),
            strip.placement   = "outside",
            legend.text       = element_text(size = textsize),
            legend.title      = element_text(size = textsize),
            aspect.ratio      = facet_ratio,
            axis.ticks        = element_line(color = "black"),
            # nature_style toggle: rectangular border vs. axis lines only
            panel.border      = if (isTRUE(nature_style)) element_blank()
                                else element_rect(color = "black", fill = NA, linewidth = 0.5),
            panel.background  = if (isTRUE(nature_style)) element_blank()
                                else element_rect(fill = NA),
            axis.line         = if (isTRUE(nature_style))
                                  element_line(color = "black", linewidth = 0.4)
                                else element_blank(),
            panel.spacing     = unit(panel_spacing, "lines"),
            axis.text         = element_text(color = "black", size = textsize),
            axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.title        = element_text(size = textsize)
        )

    dodge_width <- 0.75

    # ---- precompute per-donor connect-line geometry (non-dodged path only) ---
    # Bake .x_pos into df BEFORE the layers are added, so both geom_line and
    # geom_point read the same x. See ledger C7 for the rationale (never use
    # position_dodge on geom_line).
    use_donor_lines <- FALSE
    if (!is.null(connect_donors)) {
        if (isTRUE(dodge) && !is.null(point_color_by)) {
            warning("plotPbActExprsDiff: connect_donors is ignored in dodged ",
                    "mode (dodge = TRUE + point_color_by). Lines would not ",
                    "align with dodged point positions.", call. = FALSE)
        } else {
            use_donor_lines <- TRUE
            # Integer factor position of each row's act_cluster (matches
            # ggplot's discrete x scale: 1, 2, 3 ...)
            x_factor      <- factor(df$act_cluster,
                                    levels = levels(df$act_cluster))
            df$.x_base    <- as.numeric(x_factor)

            # Deterministic jitter (shared by lines and points)
            if (isTRUE(connect_jitter_width > 0)) {
                old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
                                get(".Random.seed", envir = .GlobalEnv) else NULL
                set.seed(connect_jitter_seed)
                df$.x_jitter <- runif(nrow(df),
                                      -connect_jitter_width,
                                       connect_jitter_width)
                if (!is.null(old_seed))
                    assign(".Random.seed", old_seed, envir = .GlobalEnv)
            } else {
                df$.x_jitter <- 0
            }
            df$.x_pos <- df$.x_base + df$.x_jitter

            # Per-donor group key: one line per (pheno_cluster × antigen × donor)
            df$.donor_group <- paste0(
                as.character(df$pheno_cluster), "\037",
                as.character(df$antigen),       "\037",
                as.character(df[[connect_donors]])
            )

            # Re-feed updated df to the ggplot via %+%
            p <- p %+% df
        }
    }

    if (dodge && !is.null(point_color_by)) {
        # ---- dodged mode: boxes and points both split by point_color_by --------------
        # Both layers share one fill scale (point_color_by). User adds one
        # scale_fill_manual(values = point_pal).
        if (geom == "boxes") {
            p <- p + geom_boxplot(
                aes(fill = .data[[point_color_by]]),
                position = position_dodge(width = dodge_width),
                color = "black", width = 0.6, linewidth = 0.3,
                alpha = 0.9, outlier.color = NA
            )
        } else {
            p <- p + stat_summary(
                fun = fun, geom = "bar",
                aes(fill = .data[[point_color_by]]),
                color = "black", position = position_dodge(width = dodge_width), alpha = 1
            )
        }
        p <- p + geom_quasirandom(
            aes(fill = .data[[point_color_by]]),
            dodge.width = dodge_width,
            width = 0.05, size = point_size, shape = 21,
            color = "black", stroke = 0.5, alpha = 0.8
        )
    } else {
        # ---- non-dodged mode: boxes filled by act_cluster, points optionally by point_color_by
        # When cluster_label_below = TRUE, widen boxes so adjacent
        # act_cluster boxes within a panel almost abut (matches
        # plotAimPaired). Standard 0.75 otherwise.
        .box_w <- if (isTRUE(cluster_label_below)) box_width_grouped else 0.75
        if (geom == "boxes") {
            p <- p + geom_boxplot(
                aes(fill = act_cluster),
                color = "black", width = .box_w, linewidth = 0.3,
                alpha = 0.9, outlier.color = NA
            )
        } else {
            p <- p + stat_summary(
                fun = fun, geom = "bar",
                aes(fill = act_cluster),
                color = "black", position = position_dodge(), alpha = 1
            )
        }
        # Donor lines BEFORE points so points render on top.
        if (isTRUE(use_donor_lines)) {
            p <- p + geom_line(
                aes(x = .x_pos, y = value, group = .donor_group),
                inherit.aes = FALSE,
                color       = connect_color,
                alpha       = connect_alpha,
                linewidth   = connect_linewidth,
                show.legend = FALSE
            )
        }
        # Points: use ggnewscale to keep point fill scale independent of box fill scale.
        # User adds: scale_fill_manual(values = box_pal) +   # first: boxes (act_cluster)
        #            scale_fill_manual(values = point_pal)   # second: points (point_color_by)
        if (!is.null(point_color_by)) {
            p <- p + ggnewscale::new_scale_fill()
            if (isTRUE(use_donor_lines)) {
                # Custom-position path: plain geom_point at .x_pos so points
                # share the line endpoints exactly. geom_quasirandom would
                # re-jitter them off the lines.
                p <- p + geom_point(
                    aes(x = .x_pos, y = value, fill = .data[[point_color_by]]),
                    inherit.aes = FALSE,
                    shape = 21, color = "black", stroke = 0.5,
                    size = point_size, alpha = 0.8
                )
            } else {
                p <- p + geom_quasirandom(
                    aes(fill = .data[[point_color_by]]),
                    width = 0.2, size = point_size, shape = 21,
                    color = "black", stroke = 0.5, alpha = 0.8
                )
            }
        } else {
            if (isTRUE(use_donor_lines)) {
                p <- p + geom_point(
                    aes(x = .x_pos, y = value),
                    inherit.aes = FALSE,
                    fill = "grey84", color = "black", stroke = 0.5,
                    shape = 21, size = point_size, alpha = 0.9
                )
            } else {
                p <- p + geom_quasirandom(
                    width = 0.2, size = point_size, shape = 21,
                    fill = "grey84", color = "black", stroke = 0.5, alpha = 0.8
                )
            }
        }
    }

    if (hide_x_labels) {
        p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
    }

    if (show_stats) {
        p <- p + ggpubr::stat_pvalue_manual(
            stat.test, label = "p.adj.signif",
            tip.length = 0.025, hide.ns = hide_ns,
            size = stat_size, bracket.nudge.y = nudge
        )
    }

    # Default layout: antigen on rows (left strip = y-axis label), pheno_cluster on cols.
    # swap = TRUE gives the transposed layout (pheno_cluster on rows, antigen on cols).
    # Antigen strip label = just the antigen name (no appended "<fun>
    # expression" line). y_label is still computed above in case it's
    # wanted elsewhere, but the per-row strip text stays compact.
    antigen_labeller <- labeller(antigen = function(x) x)

    # ---- cluster_label_below tweaks --------------------------------------
    # When TRUE and !swap, we suppress the column strip text (no bottom
    # strip), and instead draw the pheno_cluster labels as panel-scoped
    # geom_text below each panel with the rightmost letter anchored at
    # the centre x-axis tick (one tick per panel, drawn via the discrete
    # scale's `breaks`). axes = "x" makes the x-axis line appear on every
    # antigen row.
    .use_cluster_label <- isTRUE(cluster_label_below) && !isTRUE(swap)
    .eff_axes   <- if (.use_cluster_label) "x" else axes
    # CRITICAL: when cluster_label_below = TRUE we need switch = "both"
    # (antigen strip stays on the LEFT; pheno_cluster strip moves to the
    # BOTTOM). The bottom strip is KEPT (invisible via color = NA in the
    # theme block) PURELY to reserve vertical gtable space below each
    # panel -- that's the only region coord_cartesian(clip = "off") can
    # extend into. The actual visible labels are drawn by the geom_text
    # below; without the reserved strip row, the labels render outside
    # the plot region and get clipped (this is the same plotAimPaired
    # trick).
    .eff_switch <- if (.use_cluster_label) "both" else "y"

    if (label_parse) {
        if (!swap) {
            p <- p + ggh4x::facet_grid2(antigen ~ pheno_cluster, scales = scales, axes = .eff_axes,
                                         switch = .eff_switch, labeller = label_parsed)
        } else {
            p <- p + ggh4x::facet_grid2(pheno_cluster ~ antigen, scales = scales, axes = axes,
                                         labeller = label_parsed)
        }
    } else {
        if (!swap) {
            p <- p + ggh4x::facet_grid2(antigen ~ pheno_cluster, scales = scales, axes = .eff_axes,
                                         switch = .eff_switch,
                                         labeller = antigen_labeller)
        } else {
            p <- p + ggh4x::facet_grid2(pheno_cluster ~ antigen, scales = scales, axes = axes,
                                         independent = "y")
        }
    }

    # ---- cluster_label_below: single centre tick + rotated labels -------
    if (.use_cluster_label) {
        # Centre x-axis tick: use scale_x_discrete with breaks = middle
        # act_cluster level so only ONE tick is drawn per panel, at the
        # centre. (For even n, fall back to ggplot's default which is all
        # breaks -- and you should reconsider the layout.)
        .act_lvls <- levels(df$act_cluster)
        .n_act    <- length(.act_lvls)
        .x_breaks <- if (.n_act %% 2L == 1L) .act_lvls[(.n_act + 1L) / 2L]
                     else                     ggplot2::waiver()
        p <- p + ggplot2::scale_x_discrete(limits = .act_lvls, drop = FALSE,
                                           breaks = .x_breaks)

        # ONE geom_text label per pheno_cluster, tagged with the BOTTOM
        # antigen row only -- so labels render ONLY below the last row
        # of panels (into the bottom-strip gtable space reserved by the
        # invisible strip.text.x.bottom in the theme block below). All
        # rows above the bottom one get no labels.
        #
        # `levels()` returns NULL when the column is a character vector
        # (df$antigen is character at this point -- set from
        # rownames(expr_mat)); fall back to unique() so we don't end up
        # with a 0-row data.frame.
        .center_x <- (.n_act + 1L) / 2
        .pheno_lvls   <- if (is.factor(df$pheno_cluster)) levels(df$pheno_cluster)
                         else                              unique(df$pheno_cluster)
        .antigen_lvls <- if (is.factor(df$antigen))       levels(df$antigen)
                         else                              unique(df$antigen)
        .bottom_antigen <- tail(.antigen_lvls, 1L)
        .label_df <- data.frame(
            pheno_cluster = factor(.pheno_lvls, levels = .pheno_lvls),
            antigen       = factor(rep(.bottom_antigen, length(.pheno_lvls)),
                                    levels = .antigen_lvls),
            # y = -Inf renders at the panel BOTTOM (linear scale). With
            # coord_cartesian(clip = "off") below and vjust = 1 (anchor
            # at TOP of unrotated bbox), the rotated label extends
            # DOWN-LEFT into the bottom-strip area reserved below.
            x     = .center_x,
            y     = -Inf,
            label = .pheno_lvls,
            stringsAsFactors = FALSE
        )

        .label_size_pt <- if (is.null(cluster_label_size)) textsize else cluster_label_size

        p <- p +
            ggplot2::geom_text(
                data        = .label_df,
                mapping     = ggplot2::aes(x = x, y = y, label = label),
                inherit.aes = FALSE,
                angle       = cluster_label_angle,
                hjust       = cluster_label_hjust,
                vjust       = cluster_label_vjust,
                size        = .label_size_pt / 2.845
            ) +
            ggplot2::coord_cartesian(clip = "off") +
            ggplot2::theme(
                # Hide per-panel act-cluster tick LABELS but keep ticks
                # (so the centre tick set via breaks above is drawn).
                axis.text.x       = ggplot2::element_blank(),
                axis.title.x      = ggplot2::element_blank(),
                axis.ticks.x      = ggplot2::element_line(color = "black"),
                # Horizontal gap between adjacent pheno_cluster panels
                # (controlled by `cluster_panel_spacing_pt` param). > 0
                # introduces a visible gap in the per-panel x-axis line
                # between panels -- we accept that as the trade-off for
                # visible separation between the cluster "groups".
                panel.spacing.x   = grid::unit(cluster_panel_spacing_pt, "pt"),
                # Vertical gap between antigen panel rows. When the user
                # supplies `cluster_panel_spacing_y_pt`, override the
                # generic panel.spacing (set in lines earlier in this
                # theme block) with an explicit pt value. Falls back to
                # the inherited setting when NULL so existing callers
                # don't change behaviour.
                panel.spacing.y   = if (is.null(cluster_panel_spacing_y_pt))
                                        grid::unit(panel_spacing, "lines")
                                    else
                                        grid::unit(cluster_panel_spacing_y_pt, "pt"),
                # CRITICAL (mirrors plotAimPaired): the bottom strip text
                # stays present but INVISIBLE (color = NA). This reserves
                # the strip-row vertical space in the gtable so the
                # geom_text labels (drawn in panel coords, extending
                # downward past the panel via coord_cartesian(clip =
                # "off") + vjust = 1 + angle = 45) have somewhere inside
                # the plot region to render. element_blank() would
                # collapse the row -> labels would be clipped by the
                # plot-region boundary regardless of clip = "off" or
                # plot.margin bumps (plot.margin is OUTSIDE the plot
                # region, not inside it).
                strip.background.x  = ggplot2::element_rect(fill = NA, color = NA),
                strip.text.x.bottom = ggplot2::element_text(
                    size   = .label_size_pt,
                    angle  = cluster_label_angle,
                    hjust  = cluster_label_hjust,
                    vjust  = cluster_label_vjust,
                    color  = NA,
                    margin = ggplot2::margin(t = 4, b = 2)
                )
            )
    }

    # ---- tight_top_axis: per-facet y-axis termination at next pretty tick --
    # 2D facet (antigen × pheno_cluster). With scales = "free_y", each cell has
    # its own y-range; we override it per-cell so the axis ends just past the
    # data's topmost pretty break (plus a tiny overhang). When stats are shown,
    # we also accommodate the highest bracket y so brackets don't get clipped.
    # facetted_pos_scales applies the same per-facet scale list regardless of
    # facet orientation.
    if (isTRUE(tight_top_axis) && identical(scales, "free_y")) {
        # data max per (antigen, pheno_cluster)
        df_max <- df %>%
            dplyr::group_by(antigen, pheno_cluster) %>%
            dplyr::summarise(data_max = max(value, na.rm = TRUE), .groups = "drop")
        # stats max per (antigen, pheno_cluster), if present
        if (exists("stat.test", inherits = FALSE) &&
            !is.null(stat.test) && nrow(stat.test) > 0 &&
            all(c("antigen", "pheno_cluster", "y.position") %in% names(stat.test))) {
            stat_max <- stat.test %>%
                dplyr::group_by(antigen, pheno_cluster) %>%
                dplyr::summarise(stat_max = suppressWarnings(max(y.position, na.rm = TRUE)),
                                 .groups = "drop")
            df_max <- dplyr::left_join(df_max, stat_max,
                                       by = c("antigen", "pheno_cluster"))
        } else {
            df_max$stat_max <- NA_real_
        }
        # build per-facet scales (one per ROW of df_max).
        # facetted_pos_scales keys by the panel index (row-major of facets);
        # the simplest portable approach is to construct an unnamed list
        # whose order matches the facet order ggplot will use. ggh4x supports
        # both forms; we pass an unnamed list since multi-dim facet naming is
        # fiddly across ggh4x versions.
        # Order facets the same way ggplot does: row-major by row var then col var.
        if (!swap) {
            df_max <- df_max %>% dplyr::arrange(antigen, pheno_cluster)
        } else {
            df_max <- df_max %>% dplyr::arrange(pheno_cluster, antigen)
        }
        per_facet_scales <- lapply(seq_len(nrow(df_max)), function(i) {
            data_max <- df_max$data_max[i]
            stat_max <- df_max$stat_max[i]
            data_max <- ifelse(is.finite(data_max) && data_max > 0, data_max, 1)
            breaks_vec <- pretty(c(0, data_max * 1.001), n = 5)
            top_break  <- max(breaks_vec)
            candidate_top <- if (is.finite(stat_max)) max(top_break, stat_max) else top_break
            axis_top   <- candidate_top * (1 + tight_top_axis_overhang)
            ggplot2::scale_y_continuous(
                breaks = breaks_vec,
                limits = c(0, axis_top),
                expand = ggplot2::expansion(mult = c(0.02, 0))
            )
        })
        p <- p + ggh4x::facetted_pos_scales(y = per_facet_scales)
    }

    attr(p, "n_per_condition") <- n_per_cond_df
    attr(p, "source_data")     <- df

    # ---- cluster_separator: stash config + attach S3 class -----------------
    # The actual gtable overlay happens in print.plotPbActExprsDiff (below)
    # so that downstream `+ scale_fill_manual(...) + ...` operations still
    # work on the returned ggplot. ggplot2's `+.gg` returns the modified
    # plot without altering its class, so the "plotPbActExprsDiff" class
    # survives the chunk's modifier chain and the print method fires when
    # knitr auto-prints the final `p`.
    if (isTRUE(cluster_separator)) {
        attr(p, ".cluster_separator_cfg") <- list(
            linetype = cluster_separator_linetype,
            color    = cluster_separator_color,
            lwd      = cluster_separator_lwd
        )
        class(p) <- c("plotPbActExprsDiff", class(p))
    }

    return(p)
}


# ---------------------------------------------------------------------------
# S3 print method: applies cluster_separator overlay at render time.
#
# Only fires when the returned plot has class "plotPbActExprsDiff" -- which
# only happens when cluster_separator = TRUE was passed to the function (the
# class is conditionally attached just above). For all other callers the
# class isn't present, so the standard print.ggplot path runs.
#
# Flow:
#   1. ggplotGrob(x) builds the standard gtable from the (modified) ggplot
#      -- so any user-added scales / titles / ylim from the chunk's `+` chain
#      are already baked in.
#   2. We find the panel layout cells (named "panel-{row}-{col}") and identify
#      the gap layout-columns between adjacent panel columns.
#   3. For each gap column, we add a grid::segmentsGrob spanning rows
#      [min(panel.t), max(panel.b)] -- i.e. from the top of the topmost panel
#      row down through every panel.spacing.y row to the bottom of the
#      bottommost panel row. That gives a SINGLE continuous dashed line,
#      uninterrupted by the inter-row gutters.
#   4. grid.newpage() + grid.draw(g) renders the modified gtable to the
#      active device (knitr captures whatever lands on the device).
#
# We bypass patchwork::wrap_elements() entirely here: an earlier attempt
# wrapped the modified gtable via wrap_elements, but patchwork's compositor
# appeared to drop the added cluster_separator_* grobs when rebuilding its
# own outer gtable for rendering. Drawing the gtable directly avoids that.
print.plotPbActExprsDiff <- function(x, newpage = is.null(vp), vp = NULL, ...) {
    cfg <- attr(x, ".cluster_separator_cfg")
    if (is.null(cfg)) {
        # Defensive: class was attached but config missing -- fall through.
        NextMethod()
        return(invisible(x))
    }

    # IMPORTANT: grid.newpage() must come BEFORE ggplotGrob().
    # ggplotGrob -> ggplot_build -> grid-based text measurements that
    # implicitly initialize the active device's first page. If we then
    # call grid.newpage() afterwards, on cairo_pdf / pdf devices it
    # advances past that initialized page and we end up with a blank
    # page 1 + the actual plot on page 2 (this is exactly what bit f2()
    # output). Reordering to match print.ggplot (newpage -> build ->
    # draw) keeps everything on one page. Signature/vp handling also
    # mirrored on print.ggplot so callers that pass `newpage = FALSE`
    # or a `vp = ...` (e.g. when embedding into another layout) work.
    if (newpage) grid::grid.newpage()

    g <- ggplot2::ggplotGrob(x)
    panel_layout <- g$layout[grepl("^panel-", g$layout$name), , drop = FALSE]

    if (nrow(panel_layout) > 0L && length(unique(panel_layout$l)) >= 2L) {
        panel_cols <- sort(unique(panel_layout$l))
        top_row    <- min(panel_layout$t)
        bot_row    <- max(panel_layout$b)

        for (i in seq_len(length(panel_cols) - 1L)) {
            # Between adjacent panel columns in ggplot2's gtable there are
            # typically THREE layout columns: axis-r-*-i (right-axis chrome
            # of panel i), the panel.spacing.x slot, and axis-l-*-(i+1)
            # (left-axis chrome of panel i+1). The previous heuristic of
            # `panel_cols[i] + 1` landed on axis-r, which is zero-width
            # and clipped when `axes = "x"` suppresses right-axis content
            # -> the line never rendered. The integer midpoint between
            # adjacent panel cols hits the panel.spacing.x slot (e.g. for
            # panel cols 8, 12 the midpoint is 10), which is where the
            # visible gap actually lives.
            gap_col <- (panel_cols[i] + panel_cols[i + 1L]) %/% 2L
            line_grob <- grid::segmentsGrob(
                x0 = grid::unit(0.5, "npc"), y0 = grid::unit(0, "npc"),
                x1 = grid::unit(0.5, "npc"), y1 = grid::unit(1, "npc"),
                gp = grid::gpar(
                    col = cfg$color,
                    lty = cfg$linetype,
                    lwd = cfg$lwd
                )
            )
            g <- gtable::gtable_add_grob(
                g, line_grob,
                t = top_row, b = bot_row,
                l = gap_col, r = gap_col,
                z = Inf,
                name = paste0("cluster_separator_", i)
            )
        }
    }

    if (is.null(vp)) {
        grid::grid.draw(g)
    } else {
        if (is.character(vp)) grid::seekViewport(vp)
        else                  grid::pushViewport(vp)
        grid::grid.draw(g)
        grid::upViewport()
    }

    invisible(x)
}
