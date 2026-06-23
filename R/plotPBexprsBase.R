# plotPBexprsBase.R
#
# Base-R strip-plot for pseudobulk marker expression — same aesthetics as
# plotAbundancesBase.R, but the y-axis shows per-cluster median (or mean)
# expression of a single marker rather than cluster proportions.
#
# Data pipeline:
#   • For each (cluster_id × sample_id) pair, aggregate single-cell expression
#     values using fun ("median" or "mean").
#   • Metadata join, condition construction, eTreg handling, and split_col logic
#     are identical to plotAbundancesBase.
#   • Layout, jitter, segments, legend are identical to plotAbundancesBase.
#
# ─────────────────────────────────────────────────────────────────────────────
# TYPICAL CALL
# ─────────────────────────────────────────────────────────────────────────────
#
#   plotPBexprsBase(
#     x                 = sce_plot,
#     marker            = "IFNg",
#     k                 = "merging1",
#     group_by          = "T cell",
#     group_order       = c("Naive", "EM Th1", "CTL", "cTfh"),
#     split_col         = "B cell",
#     shape_by          = "patient_id",
#     meta_df           = md,
#     dodge_w           = 20,
#     group_gap         = 3,
#     per_panel_width   = 4,
#     per_panel_height  = 3.5,
#     cmv_gap           = 1.5,
#     jitter_w          = 1,
#     patient_shapes    = my_patient_shapes2,
#     stim_col          = "stim",
#     cmv_col           = "cmv",
#     n_cols            = 4L,
#     point_size        = 50,
#     point_alpha       = 1,
#     line_alpha        = 0.35,
#     line_width        = 0.55,
#     base_font_size    = 20,
#     etreg_col         = "Treg",
#     etreg_val         = "eTreg",
#     etreg_no_treg_val = "_",
#     etreg_donors      = c("LP1", "LP9", "LP10"),
#     output_file       = "AIM_TB_IFNg_expression.pdf"
#   )
#
# ─────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(SingleCellExperiment)
library(SummarizedExperiment)

# ── Visual constants (identical to plotAbundancesBase) ────────────────────────

COND_COLORS_PBE <- c(
  "DMSO CMV-" = "#999999",
  "pp65 CMV-" = "#4477BB",
  "DMSO CMV+" = "#999999",
  "pp65 CMV+" = "#CC4444",
  "eTreg"     = "#AA44BB",
  "ctfh"      = "#44AA88"
)

COND_ORDER_PBE <- c("DMSO CMV-", "pp65 CMV-", "DMSO CMV+", "pp65 CMV+", "eTreg", "ctfh")

SEGMENT_PAIRS_PBE <- list(
  c("DMSO CMV-", "pp65 CMV-"),
  c("DMSO CMV+", "pp65 CMV+"),
  c("pp65 CMV+", "eTreg")
)

SHAPE_CYCLE_PBE <- c(21L, 22L, 23L, 24L, 25L, 21L, 22L, 23L)

REP_BORDER_COLOR_PBE <- "#333333"
REP_BORDER_LWD_PBE   <- 1.5


# ══════════════════════════════════════════════════════════════════════════════
# Internal layout helpers  (identical maths to plotAbundancesBase)
# ══════════════════════════════════════════════════════════════════════════════

.computeLayout_PBE <- function(x_labels, group_gap, dodge_w, cond_levels,
                                strip_margin = 0.1, cmv_gap = 0) {
  n_cond  <- length(cond_levels)
  spacing <- dodge_w / n_cond

  raw_pos <- (seq_len(n_cond) - 1L) * spacing
  if (cmv_gap > 0 && n_cond > 2L)
    raw_pos[-(1:2)] <- raw_pos[-(1:2)] + cmv_gap
  raw_offsets <- raw_pos - mean(raw_pos)
  pd_offsets  <- setNames(raw_offsets, cond_levels)

  strip_half_w   <- max(dodge_w / 2, max(abs(raw_offsets))) + strip_margin
  centre_spacing <- 2 * strip_half_w + group_gap

  centers <- setNames(
    strip_half_w + (seq_along(x_labels) - 1L) * centre_spacing,
    x_labels
  )

  layout_df <- data.frame(
    x_label     = x_labels,
    x_centre    = unname(centers[x_labels]),
    strip_xmin  = unname(centers[x_labels]) - strip_half_w,
    strip_xmax  = unname(centers[x_labels]) + strip_half_w,
    strip_index = seq_along(x_labels),
    stringsAsFactors = FALSE
  )

  list(centers = centers, pd_offsets = pd_offsets, layout_df = layout_df)
}

.addXPos_PBE <- function(df, centers, pd_offsets, jitter_w = 0, seed = 42L) {
  df$x_centre <- unname(centers[as.character(df$x_label)])
  df$x_pos    <- df$x_centre + unname(pd_offsets[as.character(df$condition)])
  if (jitter_w > 0) {
    set.seed(seed)
    df$x_pos <- df$x_pos + runif(nrow(df), -jitter_w / 2, jitter_w / 2)
  }
  df
}

.buildSegments_PBE <- function(pts) {
  join_keys <- c("cluster_id", "patient_id", "x_label")
  if ("exp"       %in% names(pts)) join_keys <- c(join_keys, "exp")
  if ("split_val" %in% names(pts)) join_keys <- c(join_keys, "split_val")
  if ("is_rep"    %in% names(pts)) join_keys <- c(join_keys, "is_rep")

  rows <- lapply(SEGMENT_PAIRS_PBE, function(pair) {
    cf <- pair[1]; ct <- pair[2]
    if (!cf %in% pts$condition || !ct %in% pts$condition) return(NULL)

    from_df <- pts[pts$condition == cf, c(join_keys, "x_pos", "Freq"), drop = FALSE]
    to_df   <- pts[pts$condition == ct, c(join_keys, "x_pos", "Freq"), drop = FALSE]
    names(from_df)[match(c("x_pos", "Freq"), names(from_df))] <- c("seg_x0", "seg_y0")
    names(to_df  )[match(c("x_pos", "Freq"), names(to_df  ))] <- c("seg_x1", "seg_y1")

    m <- merge(from_df, to_df, by = join_keys)
    if (nrow(m) == 0L) NULL else m
  })

  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal drawing helpers
# ══════════════════════════════════════════════════════════════════════════════

.drawFacet_PBE <- function(c_pts, c_segs, layout_df, cluster,
                             shape_map, x_lim,
                             point_size, point_alpha, line_alpha, line_width,
                             font_cex, pd_offsets, cond_levels,
                             marker, fun) {

  y_vals <- if (nrow(c_pts) > 0L) c_pts$Freq else numeric(0)
  y_max  <- if (length(y_vals) > 0L && any(is.finite(y_vals))) max(y_vals, na.rm = TRUE) else 1.0
  y_min  <- if (length(y_vals) > 0L && any(is.finite(y_vals))) min(y_vals, na.rm = TRUE) else 0.0
  if (!is.finite(y_max)) y_max <- 1.0
  if (!is.finite(y_min)) y_min <- 0.0
  y_range <- max(y_max - y_min, 1e-6)
  y_pad   <- y_range * 0.08
  y_lim   <- c(y_min - y_pad, y_max + y_range * 0.18 + 1e-6)

  plot(NA, xlim = x_lim, ylim = y_lim,
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = "", bty = "n")

  # Background strips
  for (i in seq_len(nrow(layout_df))) {
    bg <- if (layout_df$strip_index[i] %% 2L == 1L) "#fafafa" else "#f3f3f3"
    rect(layout_df$strip_xmin[i], y_lim[1],
         layout_df$strip_xmax[i], y_lim[2],
         col = bg, border = NA, xpd = FALSE)
  }

  # Connecting segments
  if (!is.null(c_segs) && nrow(c_segs) > 0L) {
    has_rep_col <- "is_rep" %in% names(c_segs)

    main_segs <- if (has_rep_col) c_segs[!c_segs$is_rep, , drop = FALSE] else c_segs
    rep_segs  <- if (has_rep_col) c_segs[ c_segs$is_rep, , drop = FALSE] else data.frame()

    if (nrow(main_segs) > 0L)
      segments(main_segs$seg_x0, main_segs$seg_y0,
               main_segs$seg_x1, main_segs$seg_y1,
               col = adjustcolor("black", alpha.f = line_alpha),
               lwd = line_width * 2.5, lty = 1L)

    if (nrow(rep_segs) > 0L)
      segments(rep_segs$seg_x0, rep_segs$seg_y0,
               rep_segs$seg_x1, rep_segs$seg_y1,
               col = adjustcolor(REP_BORDER_COLOR_PBE, alpha.f = line_alpha),
               lwd = line_width * 2.0, lty = 2L)
  }

  # Points
  if (nrow(c_pts) > 0L) {
    pt_cex <- sqrt(point_size / 20)
    for (i in seq_len(nrow(c_pts))) {
      r      <- c_pts[i, ]
      fill   <- COND_COLORS_PBE[as.character(r$condition)]
      if (is.na(fill)) fill <- "#888888"
      pch    <- shape_map[[as.character(r$patient_id)]]
      if (is.null(pch) || is.na(pch)) pch <- 21L
      is_rep      <- isTRUE(r$is_rep)
      is_fillable <- pch %in% 21:25
      border <- if (is_fillable) {
        if (is_rep) REP_BORDER_COLOR_PBE else "white"
      } else {
        fill
      }
      lwd <- if (is_fillable) {
        if (is_rep) REP_BORDER_LWD_PBE else 0.5
      } else 1.5
      points(r$x_pos, r$Freq,
             pch = pch, cex = pt_cex,
             col = adjustcolor(border, alpha.f = point_alpha),
             bg  = adjustcolor(if (is_fillable) fill else NA, alpha.f = point_alpha),
             lwd = lwd)
    }
  }

  axis(2, cex.axis = font_cex * 0.8, las = 1, lwd = 0.7)
  # Y-axis label: marker name + aggregation function
  y_lab <- paste0(marker, "  [", fun, " exprs]")
  mtext(y_lab, side = 2, line = 2.8, cex = font_cex * 0.75)

  cond_x   <- unlist(lapply(layout_df$x_centre,
                            function(xc) xc + unname(pd_offsets[cond_levels])))
  cond_lab <- rep(cond_levels, nrow(layout_df))

  axis(1, at = cond_x,              labels = FALSE, tck = -0.015, lwd = 0.5)
  axis(1, at = layout_df$x_centre, labels = FALSE, tck = -0.030, lwd = 0.7)

  usr    <- par("usr")
  plot_h <- diff(usr[3:4])

  text(cond_x, usr[3] - plot_h * 0.08,
       labels = cond_lab,
       srt = 45, adj = 1, xpd = NA, cex = font_cex * 0.65)

  text(layout_df$x_centre, usr[3] - plot_h * 0.46,
       labels = layout_df$x_label,
       srt = 0, adj = 0.5, xpd = NA, cex = font_cex * 0.85, font = 2L)

  box(bty = "l", lwd = 0.7)
  title(main = cluster, cex.main = font_cex, font.main = 2, line = 0.4)
}


.drawLegend_PBE <- function(present_conds, present_patients, shape_map, font_cex,
                              has_rep = FALSE) {
  old_mar <- par(mar = c(1, 1, 1, 1))
  on.exit(par(mar = old_mar))
  plot.new()

  nc <- length(present_conds)
  np <- length(present_patients)
  sorted_pats <- sort(present_patients)

  leg_pch <- c(
    rep(22L, nc),
    vapply(sorted_pats, function(p) {
      v <- shape_map[[p]]
      if (is.null(v) || is.na(v)) 21L else as.integer(v)
    }, integer(1L))
  )

  leg_col   <- c(rep("white", nc), rep("gray40", np))
  leg_pt.bg <- c(unname(COND_COLORS_PBE[present_conds]), rep(NA_character_, np))
  leg_lwd   <- c(rep(0.5, nc), rep(NA, np))
  leg_lab   <- c(present_conds, sorted_pats)

  if (has_rep) {
    leg_lab   <- c(leg_lab,   "— Main text data",    "-- Replicate data")
    leg_pch   <- c(leg_pch,   22L,                    22L)
    leg_col   <- c(leg_col,   "white",                REP_BORDER_COLOR_PBE)
    leg_pt.bg <- c(leg_pt.bg, "gray60",               "gray60")
    leg_lwd   <- c(leg_lwd,   0.5,                    REP_BORDER_LWD_PBE)
  }

  legend("center",
    legend = leg_lab,
    pch    = leg_pch,
    col    = leg_col,
    pt.bg  = leg_pt.bg,
    pt.lwd = leg_lwd,
    pt.cex = font_cex * 1.5,
    cex    = font_cex * 0.95,
    bty    = "n"
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal: SCE → pseudobulk expression data frame
# ══════════════════════════════════════════════════════════════════════════════

.processSCE_PBE <- function(
  x,
  k,
  marker,
  assay,
  fun,
  group_by,
  shape_by,
  split_col,
  stim_col,
  cmv_col,
  etreg_col,
  etreg_val,
  etreg_donors,
  meta_df,
  use_etreg
) {
  # ── 1. Aggregate expression per (cluster_id × sample_id) ──────────────────
  c_ids <- .wl_cluster_ids(x, k)
  s_ids <- .wl_sample_ids(x)

  exprs_mat <- SummarizedExperiment::assay(x, assay)
  if (!marker %in% rownames(exprs_mat))
    stop("Marker '", marker, "' not found in rownames of assay '", assay, "'.")

  marker_exprs <- as.numeric(exprs_mat[marker, ])
  agg_fun      <- match.fun(fun)

  cell_df <- data.frame(
    cluster_id   = as.character(c_ids),
    sample_id    = as.character(s_ids),
    marker_value = marker_exprs,
    stringsAsFactors = FALSE
  )

  df <- cell_df %>%
    dplyr::group_by(cluster_id, sample_id) %>%
    dplyr::summarise(Freq = agg_fun(marker_value), .groups = "drop") %>%
    dplyr::ungroup() %>%
    as.data.frame(stringsAsFactors = FALSE)

  # ── 2. Metadata join ───────────────────────────────────────────────────────
  if (!is.null(meta_df)) {
    meta <- as.data.frame(meta_df)
    meta$sample_id <- as.character(meta$sample_id)
    vars_needed <- unique(c("sample_id", group_by, shape_by, split_col,
                             stim_col, cmv_col, "patient_id", "exp", etreg_col))
    vars_needed <- vars_needed[vars_needed %in% names(meta)]
    m_dedup <- meta %>%
      dplyr::select(dplyr::all_of(vars_needed)) %>%
      dplyr::distinct(sample_id, .keep_all = TRUE)
    df <- df %>% dplyr::left_join(m_dedup, by = "sample_id")
  } else {
    m_idx    <- match(df$sample_id, x$sample_id)
    vars_add <- unique(c(group_by, shape_by, split_col,
                         stim_col, cmv_col, "patient_id", "exp", etreg_col))
    for (v in vars_add)
      if (v %in% names(SingleCellExperiment::colData(x)))
        df[[v]] <- x[[v]][m_idx]
  }

  # ── 3. Condition variable (identical logic to plotAbundancesBase) ──────────
  if (!is.null(stim_col) && stim_col %in% names(df) &&
      !is.null(cmv_col)  && cmv_col  %in% names(df)) {

    df_raw <- df

    df <- df %>%
      dplyr::filter(
        .data[[stim_col]] %in% c("DMSO control", "pp65 AIM"),
        .data[[cmv_col]]  %in% c("neg", "pos"),
        if (use_etreg && etreg_col %in% names(df))
          .data[[etreg_col]] != etreg_val
        else
          TRUE
      ) %>%
      dplyr::mutate(
        condition = dplyr::case_when(
          .data[[stim_col]] == "DMSO control" & .data[[cmv_col]] == "neg" ~ "DMSO CMV-",
          .data[[stim_col]] == "pp65 AIM"     & .data[[cmv_col]] == "neg" ~ "pp65 CMV-",
          .data[[stim_col]] == "DMSO control" & .data[[cmv_col]] == "pos" ~ "DMSO CMV+",
          .data[[stim_col]] == "pp65 AIM"     & .data[[cmv_col]] == "pos" ~ "pp65 CMV+",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(!is.na(condition))

    if (use_etreg && etreg_col %in% names(df_raw)) {
      etreg_rows <- df_raw %>%
        dplyr::filter(
          .data[[stim_col]]  == "pp65 AIM",
          .data[[cmv_col]]   == "pos",
          .data[[etreg_col]] == etreg_val,
          patient_id         %in% etreg_donors
        ) %>%
        dplyr::mutate(condition = "eTreg")
      df <- dplyr::bind_rows(df, etreg_rows)
    }

    unmatched <- df_raw %>%
      dplyr::filter(!sample_id %in% df$sample_id) %>%
      dplyr::filter(
        !(.data[[stim_col]] %in% c("DMSO control", "pp65 AIM") &
          .data[[cmv_col]]  %in% c("neg", "pos"))
      ) %>%
      dplyr::mutate(condition = as.character(.data[[stim_col]]))
    if (nrow(unmatched) > 0L)
      df <- dplyr::bind_rows(df, unmatched)

  } else {
    stop("stim_col and cmv_col must both be present in colData / meta_df.")
  }

  # ── 4. x_label column ─────────────────────────────────────────────────────
  if (!group_by %in% names(df))
    stop("group_by column '", group_by, "' not found in data after metadata join.")
  df$x_label <- as.character(df[[group_by]])

  # ── 5. split_val column ───────────────────────────────────────────────────
  if (!is.null(split_col) && split_col %in% names(df)) {
    df$split_val <- as.character(df[[split_col]])
  } else {
    df$split_val <- "all"
  }

  # ── 6. Normalise patient_id / cluster_id ─────────────────────────────────
  if (!is.null(shape_by) && shape_by %in% names(df) && shape_by != "patient_id")
    df$patient_id <- as.character(df[[shape_by]])
  df$patient_id <- as.character(df$patient_id)
  df$cluster_id <- as.character(df$cluster_id)

  df
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal: plotting engine (expression variant of .plotFromDF)
# ══════════════════════════════════════════════════════════════════════════════

.plotFromDF_PBE <- function(
  pts_all,
  marker,
  fun,
  output_file      = NULL,
  dodge_w          = 10.0,
  group_gap        = 2.0,
  jitter_w         = 0.5,
  n_cols           = 2L,
  per_panel_width  = 11.0,
  per_panel_height = 3.5,
  base_font_size   = 9.0,
  point_size       = 35.0,
  point_alpha      = 0.70,
  line_alpha        = 0.35,
  line_width        = 0.55,
  split_val_filter = NULL,
  group_order      = NULL,
  clusters_order   = NULL,
  patient_shapes   = NULL,
  cmv_gap          = 0
) {
  pts_all$cluster_id <- as.character(pts_all$cluster_id)
  pts_all$patient_id <- as.character(pts_all$patient_id)
  if (!"is_rep" %in% names(pts_all)) pts_all$is_rep <- FALSE

  has_rep <- any(pts_all$is_rep, na.rm = TRUE)

  cond_levels    <- intersect(COND_ORDER_PBE, unique(pts_all$condition))
  all_x_labels   <- sort(unique(as.character(pts_all$x_label)))
  ordered_labels <- if (!is.null(group_order)) {
    c(intersect(group_order, all_x_labels), setdiff(all_x_labels, group_order))
  } else all_x_labels

  lay        <- .computeLayout_PBE(ordered_labels, group_gap, dodge_w, cond_levels,
                                    cmv_gap = cmv_gap)
  centers    <- lay$centers
  pd_offsets <- lay$pd_offsets
  layout_df  <- lay$layout_df

  x_margin <- 0.5
  x_lim    <- c(min(layout_df$strip_xmin) - x_margin,
                max(layout_df$strip_xmax) + x_margin)

  pts_all <- .addXPos_PBE(pts_all, centers, pd_offsets, jitter_w = jitter_w)

  all_patients <- sort(unique(pts_all$patient_id))
  if (is.null(patient_shapes)) {
    shape_map <- setNames(
      as.list(SHAPE_CYCLE_PBE[
        (seq_along(all_patients) - 1L) %% length(SHAPE_CYCLE_PBE) + 1L
      ]),
      all_patients
    )
  } else {
    shape_map <- as.list(patient_shapes)
    missing_p <- setdiff(all_patients, names(shape_map))
    for (i in seq_along(missing_p))
      shape_map[[missing_p[i]]] <-
        SHAPE_CYCLE_PBE[(i - 1L) %% length(SHAPE_CYCLE_PBE) + 1L]
  }

  all_clusters     <- sort(unique(pts_all$cluster_id))
  ordered_clusters <- if (!is.null(clusters_order)) {
    c(intersect(clusters_order, all_clusters), setdiff(all_clusters, clusters_order))
  } else all_clusters

  split_vals <- sort(unique(as.character(pts_all$split_val)))
  if (!is.null(split_val_filter))
    split_vals <- intersect(split_vals, split_val_filter)

  font_cex <- base_font_size / 12

  if (!is.null(output_file)) {
    n_cl   <- length(ordered_clusters)
    n_rows <- ceiling((n_cl + 1L) / n_cols)
    fig_w  <- per_panel_width  * n_cols + 1.5
    fig_h  <- per_panel_height * n_rows
    ext    <- tolower(tools::file_ext(output_file))
    if (ext == "pdf") {
      pdf(output_file, width = fig_w, height = fig_h)
    } else {
      png(output_file,
          width  = round(fig_w * 300),
          height = round(fig_h * 300),
          res    = 300)
    }
    on.exit(dev.off(), add = TRUE)
  }

  for (sv in split_vals) {
    sv_pts <- pts_all[as.character(pts_all$split_val) == sv, , drop = FALSE]
    if (nrow(sv_pts) == 0L) next

    present_clusters <- intersect(ordered_clusters, unique(sv_pts$cluster_id))
    if (length(present_clusters) == 0L) next

    sv_segs  <- .buildSegments_PBE(sv_pts)
    n_panels <- length(present_clusters)
    n_cells  <- n_panels + 1L
    n_rows   <- ceiling(n_cells / n_cols)

    old_par <- par(
      mfrow = c(n_rows, n_cols),
      mar   = c(9.5, 4.5, 2.0, 0.5),
      oma   = c(0, 0, if (sv != "all") 2.5 else 0, 0)
    )
    on.exit(par(old_par), add = TRUE)

    for (cluster in present_clusters) {
      c_pts  <- sv_pts[sv_pts$cluster_id == cluster, , drop = FALSE]
      c_segs <- if (nrow(sv_segs) > 0L)
                  sv_segs[sv_segs$cluster_id == cluster, , drop = FALSE]
                else data.frame()
      .drawFacet_PBE(c_pts, c_segs, layout_df, cluster,
                     shape_map, x_lim,
                     point_size, point_alpha, line_alpha, line_width, font_cex,
                     pd_offsets, cond_levels,
                     marker = marker, fun = fun)
    }

    present_conds    <- intersect(COND_ORDER_PBE, unique(sv_pts$condition))
    present_patients <- sort(unique(sv_pts$patient_id))
    .drawLegend_PBE(present_conds, present_patients, shape_map, font_cex,
                    has_rep = has_rep)

    n_empty <- n_rows * n_cols - n_cells
    if (n_empty > 0L) for (i in seq_len(n_empty)) { par(mar=c(0,0,0,0)); plot.new() }

    if (sv != "all")
      mtext(sv, outer = TRUE, cex = font_cex * 1.3, font = 2, line = 1)
  }

  invisible(NULL)
}


# ══════════════════════════════════════════════════════════════════════════════
# plotPBexprsBase  — main entry point
# ══════════════════════════════════════════════════════════════════════════════

plotPBexprsBase <- function(
  x,
  marker,                          # REQUIRED: marker name (rowname in assay)
  k                 = "meta20",
  assay             = "exprs",     # which SummarizedExperiment assay to use
  fun               = "median",    # "median" or "mean"
  group_by          = "condition",
  shape_by          = NULL,
  split_col         = NULL,
  n_cols            = 2L,
  clusters_order    = NULL,
  group_order       = NULL,
  meta_df           = NULL,
  stim_col          = "stim",
  cmv_col           = "cmv",
  dodge_w           = 10.0,
  group_gap         = 2.0,
  jitter_w          = 0.5,
  line_alpha        = 0.35,
  line_width        = 0.55,
  etreg_col         = "Treg",
  etreg_no_treg_val = "_",
  etreg_val         = "eTreg",
  etreg_donors      = NULL,
  output_file       = NULL,
  per_panel_width   = 11.0,
  per_panel_height  = 3.5,
  base_font_size    = 9.0,
  point_size        = 35.0,
  point_alpha       = 0.70,
  split_val_filter  = NULL,
  patient_shapes    = NULL,
  cmv_gap           = 0,
  # ---- f2-style auto-save (added 2026-05-25; mirrors plotAbundancesBase) ---
  # When `title` is non-NULL, auto-saves to here::here(<dated folder>) using
  # the same convention as f2() / saveFig() / plotAbundancesBase. PDF by
  # default; pass `out_format = "svg"` for SVG. `saveExcel = TRUE` writes
  # the pseudo-bulk dataframe alongside.
  title             = NULL,
  saveExcel         = FALSE,
  out_format        = c("pdf", "svg")
) {
  fun <- match.arg(fun, c("median", "mean"))
  out_format <- match.arg(out_format)

  .wl_check_sce(x, TRUE)
  k <- .wl_check_k(x, k)

  use_etreg <- !is.null(etreg_donors)

  # ---- f2-style folder + filename construction ---------------------------
  excel_filename <- NULL
  if (!is.null(title)) {
    if (!is.null(output_file)) {
      warning("plotPBexprsBase: both `title` and `output_file` provided. ",
              "`output_file` takes precedence for the figure file; ",
              "`saveExcel` will use `title` for the data filename.")
    }
    if (!requireNamespace("rmdhelp", quietly = TRUE)) {
      rmd_name_noext <- "plotPBexprsBase"
    } else {
      rmd_file_path <- tryCatch(rmdhelp::get_this_rmd_file(),
                                error = function(e) NULL)
      rmd_name_noext <- if (!is.null(rmd_file_path))
        tools::file_path_sans_ext(basename(rmd_file_path)) else "plotPBexprsBase"
    }
    folder_name <- here::here(paste0(Sys.Date(), "_", rmd_name_noext))
    if (!dir.exists(folder_name)) dir.create(folder_name, recursive = TRUE)
    timestamp <- base::format(Sys.time(), "%Y-%m-%d_%H.%M.%S")
    if (is.null(output_file)) {
      output_file <- file.path(folder_name,
                               paste0(timestamp, "_", title, ".", out_format))
    }
    if (isTRUE(saveExcel)) {
      excel_filename <- file.path(folder_name,
                                  paste0(timestamp, "_", title, "_data.xlsx"))
    }
  }

  pts_all        <- .processSCE_PBE(
    x            = x,
    k            = k,
    marker       = marker,
    assay        = assay,
    fun          = fun,
    group_by     = group_by,
    shape_by     = shape_by,
    split_col    = split_col,
    stim_col     = stim_col,
    cmv_col      = cmv_col,
    etreg_col    = etreg_col,
    etreg_val    = etreg_val,
    etreg_donors = etreg_donors,
    meta_df      = meta_df,
    use_etreg    = use_etreg
  )
  pts_all$is_rep <- FALSE

  .plotFromDF_PBE(
    pts_all          = pts_all,
    marker           = marker,
    fun              = fun,
    output_file      = output_file,
    dodge_w          = dodge_w,
    group_gap        = group_gap,
    jitter_w         = jitter_w,
    n_cols           = as.integer(n_cols),
    per_panel_width  = per_panel_width,
    per_panel_height = per_panel_height,
    base_font_size   = base_font_size,
    point_size       = point_size,
    point_alpha      = point_alpha,
    line_alpha       = line_alpha,
    line_width       = line_width,
    split_val_filter = split_val_filter,
    group_order      = group_order,
    clusters_order   = clusters_order,
    patient_shapes   = patient_shapes,
    cmv_gap          = cmv_gap
  )

  # ── Save underlying data to Excel (if requested) ───────────────────────
  if (!is.null(excel_filename)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      warning("openxlsx is not installed; skipping Excel data save.")
    } else {
      wb <- openxlsx::createWorkbook()
      openxlsx::modifyBaseFont(wb, fontSize = 12, fontColour = "black",
                               fontName = "Arial")
      openxlsx::addWorksheet(wb, "data")
      openxlsx::writeData(wb, "data", as.data.frame(pts_all))
      if (ncol(pts_all) >= 1) {
        openxlsx::setColWidths(wb, "data", cols = seq_len(ncol(pts_all)),
                               widths = "auto")
        openxlsx::freezePane(wb, "data", firstRow = TRUE)
        openxlsx::addFilter(wb, "data", rows = 1,
                            cols = seq_len(ncol(pts_all)))
      }
      openxlsx::saveWorkbook(wb, excel_filename, overwrite = TRUE)
      message("Saved data to: ", excel_filename)
    }
  }
  if (!is.null(output_file)) {
    message("Saved figure to: ", output_file)
  }

  invisible(NULL)
}

# Give plotPBexprsBase access to CATALYST internals (.check_sce, .check_k)
