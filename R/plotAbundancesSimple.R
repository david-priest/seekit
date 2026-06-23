# plotAbundancesSimple.R
#
# Simplified strip-plot for B cell abundance: no Tregs, no CMV split.
# Naive B cell T cell conditions on the left half of the x-axis;
# Memory B cell T cell conditions on the right half.
#
# Each column = one (B cell type × T cell subset) combination.
# Colors  = T cell subset (consistent across Naive and Memory halves).
# Shapes  = patient/donor (same pch 21-25 cycle as plotAbundancesBase).
# Lines   = optional; connect same patient within each B cell half only —
#           never drawn across the Naive / Memory boundary.
# Median  = crossbar drawn for every column.
#
# ─────────────────────────────────────────────────────────────────────────────
# TYPICAL CALL
# ─────────────────────────────────────────────────────────────────────────────
#
#   plotAbundancesSimple(
#     x            = sce,
#     k            = "meta20",
#     b_cell_col   = "B cell",
#     t_cell_col   = "T cell",
#     t_cell_order = c("Naive", "EM Th1", "CTL", "cTfh"),
#     draw_lines   = TRUE,
#     output_file  = "abundances_simple.pdf"
#   )
#
# ─────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(SingleCellExperiment)

# ── Default colour palette for T cell subsets ─────────────────────────────────
# (Same set of hues used in plotAbundancesBase for conditions)
.T_COLORS_SIMPLE <- c(
  "#4477BB",   # blue
  "#CC4444",   # red
  "#44AA77",   # teal
  "#AA44AA",   # purple
  "#DDAA33",   # gold
  "#77CCEE",   # sky blue
  "#EE8833",   # orange
  "#999999"    # gray
)

# Filled-symbol cycle — identical to plotAbundancesBase
# pch 21-25: col = border (white), bg = fill colour
SHAPE_CYCLE_SIMPLE <- c(21L, 22L, 23L, 24L, 25L, 21L, 22L, 23L)

# Default human-readable labels for B cell types
.B_LABELS_SIMPLE <- c(Naive = "Naive B cells", Mem = "Memory B cells")


# ══════════════════════════════════════════════════════════════════════════════
# Internal: layout helpers
# ══════════════════════════════════════════════════════════════════════════════

# Build data.frames describing x-positions of every column and every B cell
# type block (for background shading and block-header labels).
#
# Returns: list(col_df, b_block_df, x_lim)
#   col_df     : b_cell, t_cell, x_centre, col_idx, strip_half_w
#   b_block_df : b_cell, x_start, x_end, x_mid
#   x_lim      : numeric(2)
.buildSimpleLayout <- function(b_order, t_order,
                                col_gap = 1.0, b_gap = 2.5) {
  x        <- 0.0
  col_idx  <- 0L
  col_rows <- list()
  b_rows   <- list()

  for (b in b_order) {
    # Insert larger gap before every B cell group after the first
    if (col_idx > 0L) x <- x + (b_gap - col_gap)
    b_start <- x - col_gap / 2

    for (t in t_order) {
      col_idx <- col_idx + 1L
      col_rows[[col_idx]] <- data.frame(
        b_cell       = b,
        t_cell       = t,
        x_centre     = x,
        col_idx      = col_idx,
        strip_half_w = col_gap / 2,
        stringsAsFactors = FALSE
      )
      x <- x + col_gap
    }

    b_end <- x - col_gap / 2
    b_rows[[length(b_rows) + 1L]] <- data.frame(
      b_cell  = b,
      x_start = b_start,
      x_end   = b_end,
      x_mid   = (b_start + b_end) / 2,
      stringsAsFactors = FALSE
    )
  }

  col_df     <- do.call(rbind, col_rows)
  b_block_df <- do.call(rbind, b_rows)

  x_margin <- col_gap * 0.6
  x_lim    <- c(min(col_df$x_centre) - x_margin,
                max(col_df$x_centre) + x_margin)

  list(col_df = col_df, b_block_df = b_block_df, x_lim = x_lim)
}


# Add a jittered x_pos column to pts using the col_df lookup table.
.addSimpleXPos <- function(pts, col_df, jitter_w = 0.3, seed = 42L) {
  key_map <- setNames(
    col_df$x_centre,
    paste(col_df$b_cell, col_df$t_cell, sep = "|||")
  )
  pts$x_centre <- unname(key_map[paste(pts$b_cell, pts$t_cell, sep = "|||")])
  set.seed(seed)
  pts$x_pos <- pts$x_centre +
    runif(nrow(pts), -jitter_w / 2, jitter_w / 2)
  pts
}


# Build within-group patient line segments for one cluster panel.
# Lines connect adjacent T cell columns for the same patient within each
# B cell half; they never cross the Naive / Memory boundary.
.buildSimpleSegments <- function(c_pts, t_order) {
  if (nrow(c_pts) == 0L) return(data.frame())

  rows <- list()
  for (b_val in unique(c_pts$b_cell)) {
    b_pts <- c_pts[c_pts$b_cell == b_val, , drop = FALSE]
    for (pat in unique(b_pts$patient_id)) {
      p_pts <- b_pts[b_pts$patient_id == pat, , drop = FALSE]
      # Sort rows according to t_cell_order
      idx   <- match(intersect(t_order, p_pts$t_cell), p_pts$t_cell)
      p_pts <- p_pts[idx[!is.na(idx)], , drop = FALSE]
      if (nrow(p_pts) < 2L) next
      for (i in seq_len(nrow(p_pts) - 1L)) {
        rows[[length(rows) + 1L]] <- data.frame(
          x0 = p_pts$x_pos[i],     y0 = p_pts$Freq[i],
          x1 = p_pts$x_pos[i + 1L], y1 = p_pts$Freq[i + 1L],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal: drawing helpers
# ══════════════════════════════════════════════════════════════════════════════

.drawSimpleFacet <- function(c_pts, segs, col_df, b_block_df,
                              cluster, t_colors, shape_map,
                              x_lim, t_order, b_labels, draw_lines,
                              point_size, point_alpha,
                              line_alpha, line_width, font_cex,
                              y_axis_title = NULL) {

  y_max <- if (nrow(c_pts) > 0L) max(c_pts$Freq, na.rm = TRUE) else 1.0
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1.0
  y_pad <- y_max * 0.08
  y_lim <- c(-y_pad, y_max * 1.22 + 1e-4)

  plot(NA, xlim = x_lim, ylim = y_lim,
       xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = "", bty = "n")

  # Faint horizontal gridlines at major y-ticks (drawn FIRST so segments +
  # points render on top). Matches the plotAbundancesBase style — helps the
  # eye read absolute % values when n is small.
  abline(h = axTicks(2), col = "#eeeeee", lwd = 0.3)

  # Dashed vertical separator between Naive and Memory halves
  if (nrow(b_block_df) == 2L) {
    x_sep <- (b_block_df$x_end[1L] + b_block_df$x_start[2L]) / 2
    abline(v = x_sep, col = "#bbbbbb", lty = 2L, lwd = 0.8)
  }

  # Connecting lines (within B cell half, across T cell columns, same patient)
  if (draw_lines && !is.null(segs) && nrow(segs) > 0L) {
    segments(segs$x0, segs$y0, segs$x1, segs$y1,
             col = adjustcolor("black", alpha.f = line_alpha),
             lwd = line_width * 2.5)
  }


  # Data points (pch 21-25: same-colour outline + alpha'd fill, others:
  # coloured lines). Border = same hue as fill but FULL opacity; bg = same
  # hue with `point_alpha`. Gives each point a crisp coloured edge against
  # an alphad interior — replaces the previous "white halo" look that read
  # as silly when stacked on grey gridlines.
  if (nrow(c_pts) > 0L) {
    pt_cex <- sqrt(point_size / 20)
    for (i in seq_len(nrow(c_pts))) {
      r   <- c_pts[i, ]
      col <- t_colors[[as.character(r$t_cell)]]
      if (is.null(col) || is.na(col)) col <- "#888888"
      pch <- shape_map[[as.character(r$patient_id)]]
      if (is.null(pch) || is.na(pch)) pch <- 21L

      is_fillable <- pch %in% 21:25

      points(r$x_pos, r$Freq,
             pch = pch,
             cex = pt_cex,
             col = col,                                       # full-opacity edge
             bg  = adjustcolor(if (is_fillable) col else NA,
                               alpha.f = point_alpha),
             lwd = if (is_fillable) 0.7 else 1.5)             # 0.5 → 0.7 so
                                                              # the new same-
                                                              # colour border
                                                              # reads cleanly
    }
  }

  # Median crossbar for each column
  for (i in seq_len(nrow(col_df))) {
    b_val <- col_df$b_cell[i]
    t_val <- col_df$t_cell[i]
    grp   <- c_pts[c_pts$b_cell == b_val & c_pts$t_cell == t_val, , drop = FALSE]
    if (nrow(grp) == 0L) next
    med <- median(grp$Freq, na.rm = TRUE)
    xc  <- col_df$x_centre[i]
    col <- t_colors[[t_val]]
    if (is.null(col) || is.na(col)) col <- "#888888"
    segments(xc - 0.28, med, xc + 0.28, med,
             col = "black", lwd = 1.5, lend = "round")
  }

  # Y-axis. Multipliers bumped to 0.9 (was 0.8 / 0.65) so smallest text
  # remains ≥6pt after downscaling to a paper-column figure. The "Proportion
  # [%]" y-axis title is drawn ONLY when y_axis_title is non-NULL — passing
  # one title per panel is usually redundant once readers see the y-tick
  # numbers; default is off.
  axis(2, cex.axis = font_cex * 0.9, las = 1, lwd = 0.7)
  if (!is.null(y_axis_title) && nzchar(y_axis_title)) {
    mtext(y_axis_title, side = 2, line = 2.8, cex = font_cex * 0.9)
  }

  # ---- Bottom-margin labels (line-based positioning) -----------------------
  # `text()` y positions are converted to *margin-line* offsets below the
  # axis, NOT to fractions of the user-coord plot height. Fraction-based
  # positions shift with per_panel_height: when margins grow (e.g. to make
  # room for a bigger title) the plot region shrinks and the fraction-based
  # label positions collapse INTO the rotated-label projection. The fix is
  # to make the offsets constant in absolute units (lines / inches) so the
  # T-cell labels and B-cell block headers stay correctly stacked.
  usr <- par("usr")
  pin <- par("pin")             # plot region in inches (w, h)
  cin <- par("cin")             # 1 default-text line in inches (w, h)
  lh_user_y <- cin[2] * (diff(usr[3:4]) / pin[2])  # 1 line in user-y units

  # Primary x-axis: T cell subset labels (rotated 45°). Anchored 0.7 lines
  # below the axis; the label then projects down-left 5-6 lines depending
  # on label length and cex (cex font_cex × 0.9).
  axis(1, at = col_df$x_centre, labels = FALSE, tck = -0.025, lwd = 0.7)
  text(col_df$x_centre, usr[3] - 0.7 * lh_user_y,
       labels = col_df$t_cell,
       srt = 45, adj = c(1, 0.5), xpd = NA, cex = font_cex * 0.9)

  # Secondary axis: bold B cell type block headers. Placed 6 lines below
  # the axis so they sit BELOW the projection of the rotated T-cell labels
  # (which extends ~5-6 lines down-left from each tick).
  for (i in seq_len(nrow(b_block_df))) {
    lab <- b_labels[b_block_df$b_cell[i]]
    if (is.na(lab)) lab <- b_block_df$b_cell[i]
    text(b_block_df$x_mid[i], usr[3] - 6.0 * lh_user_y,
         labels = lab,
         font = 2L, xpd = NA, cex = font_cex * 0.95)
  }

  box(bty = "l", lwd = 0.7)

  # Panel title via .drawClippedMarginText (defined in plotAbundancesBase.R,
  # available in the global env when both files are sourced). Produces a
  # tight per-title PDF clip rect instead of the default panel-sized mask.
  .drawClippedMarginText(cluster, side = 3, line = 1.4,
                         cex = font_cex, font = 2L)
}


.drawSimpleLegend <- function(t_order, t_colors,
                               present_patients, shape_map, font_cex) {
  old_mar <- par(mar = c(1, 1, 1, 1))
  on.exit(par(mar = old_mar))
  plot.new()

  nt          <- length(t_order)
  sorted_pats <- sort(present_patients)
  np          <- length(sorted_pats)

  pch_vec <- c(
    rep(22L, nt),
    vapply(sorted_pats, function(p) {
      v <- shape_map[[p]]
      if (is.null(v) || is.na(v)) 21L else as.integer(v)
    }, integer(1L))
  )
  col_vec <- c(rep("white",  nt), rep("gray40", np))
  bg_vec  <- c(
    vapply(t_order, function(t) {
      co <- t_colors[[t]]
      if (is.null(co) || is.na(co)) "#888888" else co
    }, character(1L)),
    rep(NA_character_, np)
  )

  legend("center",
    legend = c(t_order, sorted_pats),
    pch    = pch_vec,
    col    = col_vec,
    pt.bg  = bg_vec,
    pt.cex = font_cex * 1.5,
    cex    = font_cex * 0.95,
    bty    = "n"
  )
}


# ══════════════════════════════════════════════════════════════════════════════
# plotAbundancesSimple  — main entry point
# ══════════════════════════════════════════════════════════════════════════════

plotAbundancesSimple <- function(
  x,                             # SingleCellExperiment
  k,                             # clustering key (e.g. "meta20")
  b_cell_col,                    # colData column for B cell type ("B cell")
  t_cell_col,                    # colData column for T cell subset ("T cell")
  meta_df           = NULL,      # optional external metadata data.frame
  patient_col       = "patient_id",
  b_cell_order      = c("Naive", "Mem"),   # left → right
  t_cell_order      = NULL,      # explicit T cell column order; NULL = alphabetical
  clusters_order    = NULL,      # ordering of B cell cluster facet panels
  t_cell_colors     = NULL,      # named character vector of hex colours; NULL = auto
  draw_lines        = FALSE,     # connect same patient within each B cell half
  n_cols            = 2L,        # panel columns in the output grid
  col_gap           = 1.0,       # centre-to-centre spacing within a B cell half
  b_gap             = 2.5,       # centre-to-centre spacing across the Naive/Mem gap
  jitter_w          = 0.3,       # total horizontal jitter width
  output_file       = NULL,      # path for PDF/PNG; NULL = draw to current device
  per_panel_width   = 4.5,       # inches per facet panel
  per_panel_height  = 3.5,
  base_font_size    = 9.0,
  point_size        = 35.0,
  point_alpha       = 0.70,
  line_alpha        = 0.35,
  patient_shapes    = NULL,             # named integer vector (patient → pch); NULL = auto
  line_width        = 0.55,
  # Per-panel y-axis title (e.g. "Proportion [%]"). NULL = no label, which
  # is the default — the y-tick numbers already convey scale and a label
  # per panel is usually visual noise. Pass a string to enable.
  y_axis_title      = NULL,
  # ---- f2-style auto-save (added 2026-05-25; mirrors plotAbundancesBase) ---
  # When `title` is non-NULL the function writes the figure into the standard
  # dated folder (here::here(<YYYY-MM-DD>_<rmd-name>/)). PDF by default;
  # pass `out_format = "svg"` for SVG. `saveExcel = TRUE` writes the
  # underlying long-format data alongside.
  #
  # Explicit `output_file = ...` still wins (backwards-compatible). If both
  # `title` and `output_file` are set, `output_file` takes precedence for
  # the figure; `saveExcel` still derives the xlsx name from `title`.
  title             = NULL,
  saveExcel         = FALSE,
  out_format        = c("pdf", "svg")
) {
  out_format <- match.arg(out_format)

  # ---- f2-style folder + filename construction ---------------------------
  excel_filename <- NULL
  if (!is.null(title)) {
    if (!is.null(output_file)) {
      warning("plotAbundancesSimple: both `title` and `output_file` provided. ",
              "`output_file` takes precedence for the figure; `saveExcel` ",
              "will use `title` for the data filename.")
    }
    if (!requireNamespace("rmdhelp", quietly = TRUE)) {
      rmd_name_noext <- "plotAbundancesSimple"
    } else {
      rmd_file_path <- tryCatch(rmdhelp::get_this_rmd_file(),
                                error = function(e) NULL)
      rmd_name_noext <- if (!is.null(rmd_file_path))
        tools::file_path_sans_ext(basename(rmd_file_path)) else "plotAbundancesSimple"
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

  # ── 1. Cluster frequencies from SCE ─────────────────────────────────────────
  .wl_check_sce(x, TRUE)
  k <- .wl_check_k(x, k)

  c_ids <- .wl_cluster_ids(x, k)
  s_ids <- .wl_sample_ids(x)
  ns    <- table(cluster_id = c_ids, sample_id = s_ids)
  fq    <- prop.table(ns, 2) * 100
  df    <- as.data.frame(fq, stringsAsFactors = FALSE) %>%
    dplyr::filter(!is.nan(Freq))

  # ── 2. Metadata join ─────────────────────────────────────────────────────────
  vars_needed <- unique(c("sample_id", b_cell_col, t_cell_col, patient_col))
  if (!is.null(meta_df)) {
    meta <- as.data.frame(meta_df)
    meta$sample_id <- as.character(meta$sample_id)
    vars_present   <- vars_needed[vars_needed %in% names(meta)]
    m_dedup <- meta %>%
      dplyr::select(dplyr::all_of(vars_present)) %>%
      dplyr::distinct(sample_id, .keep_all = TRUE)
    df <- df %>% dplyr::left_join(m_dedup, by = "sample_id")
  } else {
    m_idx <- match(df$sample_id, x$sample_id)
    for (v in setdiff(vars_needed, "sample_id"))
      if (v %in% names(SingleCellExperiment::colData(x)))
        df[[v]] <- x[[v]][m_idx]
  }

  # ── 3. Standardise column names ──────────────────────────────────────────────
  df$b_cell     <- as.character(df[[b_cell_col]])
  df$t_cell     <- as.character(df[[t_cell_col]])
  df$patient_id <- as.character(df[[patient_col]])
  df$cluster_id <- as.character(df$cluster_id)

  df <- df[!is.na(df$b_cell) & !is.na(df$t_cell), , drop = FALSE]

  # ── 4. Resolve orderings ─────────────────────────────────────────────────────
  b_cell_order <- intersect(b_cell_order, unique(df$b_cell))
  if (length(b_cell_order) == 0L)
    stop("None of b_cell_order values found in column '", b_cell_col, "'.")

  t_present <- unique(df$t_cell)
  if (is.null(t_cell_order)) {
    t_cell_order <- sort(t_present)
  } else {
    t_cell_order <- c(
      intersect(t_cell_order, t_present),
      setdiff(t_present, t_cell_order)
    )
  }

  all_clusters     <- sort(unique(df$cluster_id))
  ordered_clusters <- if (!is.null(clusters_order)) {
    c(intersect(clusters_order, all_clusters),
      setdiff(all_clusters, clusters_order))
  } else all_clusters

  # ── 5. Assign T cell colours ─────────────────────────────────────────────────
  if (is.null(t_cell_colors)) {
    t_cell_colors <- setNames(
      .T_COLORS_SIMPLE[(seq_along(t_cell_order) - 1L) %%
                         length(.T_COLORS_SIMPLE) + 1L],
      t_cell_order
    )
  } else {
    missing_t <- setdiff(t_cell_order, names(t_cell_colors))
    for (i in seq_along(missing_t))
      t_cell_colors[[missing_t[i]]] <-
        .T_COLORS_SIMPLE[(i - 1L) %% length(.T_COLORS_SIMPLE) + 1L]
  }
  t_colors_list <- as.list(t_cell_colors)

  # ── 6. Build layout ──────────────────────────────────────────────────────────
  lay        <- .buildSimpleLayout(b_cell_order, t_cell_order,
                                    col_gap = col_gap, b_gap = b_gap)
  col_df     <- lay$col_df
  b_block_df <- lay$b_block_df
  x_lim      <- lay$x_lim

  # ── 7. Jittered x positions ──────────────────────────────────────────────────
  df <- .addSimpleXPos(df, col_df, jitter_w = jitter_w)

  # ── 8. Patient → shape mapping ───────────────────────────────────────────────
  all_patients <- sort(unique(df$patient_id))
  if (is.null(patient_shapes)) {
    shape_map <- setNames(
      as.list(SHAPE_CYCLE_SIMPLE[
        (seq_along(all_patients) - 1L) %% length(SHAPE_CYCLE_SIMPLE) + 1L
      ]),
      all_patients
    )
  } else {
    shape_map <- as.list(patient_shapes)
    missing_p <- setdiff(all_patients, names(shape_map))
    for (i in seq_along(missing_p))
      shape_map[[missing_p[i]]] <-
        SHAPE_CYCLE_SIMPLE[(i - 1L) %% length(SHAPE_CYCLE_SIMPLE) + 1L]
  }

  # ── 9. B cell block labels ───────────────────────────────────────────────────
  b_labels <- .B_LABELS_SIMPLE
  for (b in b_cell_order)
    if (!b %in% names(b_labels)) b_labels[b] <- b

  font_cex <- base_font_size / 12

  # ── 10. Calculate figure dimensions ──────────────────────────────────────────
  n_cl       <- length(ordered_clusters)
  n_rows_fig <- ceiling((n_cl + 1L) / n_cols)   # +1 for legend panel
  fig_w      <- per_panel_width  * n_cols
  fig_h      <- per_panel_height * n_rows_fig

  # ── 11. Build drawing closure ─────────────────────────────────────────────────
  draw_fn <- function() {
    n_panels <- length(ordered_clusters)
    n_cells  <- n_panels + 1L
    n_rows   <- ceiling(n_cells / n_cols)

    # Margins widened to match plotAbundancesBase's geometry after the
    # multiplier bumps. bottom = T-cell labels + bold B-cell block headers
    # (drawn 32% of plot-height below the axis); top = panel title clearance
    # over the panel above's bottom labels. See margin-rationale comments
    # in plotAbundancesBase.R.
    old_par <- par(
      mfrow = c(n_rows, n_cols),
      mar   = c(8.5, 4.5, 5.5, 0.5),
      oma   = c(0, 0, 0, 0)
    )
    on.exit(par(old_par), add = TRUE)

    for (cluster in ordered_clusters) {
      c_pts <- df[df$cluster_id == cluster, , drop = FALSE]
      segs  <- if (draw_lines) .buildSimpleSegments(c_pts, t_cell_order)
               else NULL
      .drawSimpleFacet(
        c_pts        = c_pts,
        segs         = segs,
        col_df       = col_df,
        b_block_df   = b_block_df,
        cluster      = cluster,
        t_colors     = t_colors_list,
        shape_map    = shape_map,
        x_lim        = x_lim,
        t_order      = t_cell_order,
        b_labels     = b_labels,
        draw_lines   = draw_lines,
        point_size   = point_size,
        point_alpha  = point_alpha,
        line_alpha   = line_alpha,
        line_width   = line_width,
        font_cex     = font_cex,
        y_axis_title = y_axis_title
      )
    }

    .drawSimpleLegend(t_cell_order, t_colors_list,
                      all_patients, shape_map, font_cex)

    n_empty <- n_rows * n_cols - n_cells
    if (n_empty > 0L)
      for (i in seq_len(n_empty)) { par(mar = c(0, 0, 0, 0)); plot.new() }
  }

  # ── 12. Save to output_file if requested ──────────────────────────────────────
  if (!is.null(output_file)) {
    ext <- tolower(tools::file_ext(output_file))
    if (ext == "pdf") {
      pdf(output_file, width = fig_w, height = fig_h)
    } else {
      png(output_file,
          width  = round(fig_w * 300),
          height = round(fig_h * 300),
          res    = 300)
    }
    draw_fn()
    dev.off()
  }

  # ── 13. Save underlying data to Excel (if requested via title+saveExcel) ──
  # Mirrors f2()'s saveExcel = TRUE behaviour. Writes the long-format
  # per-point dataframe (one row per sample × cluster × t/b cell) to a
  # single "data" sheet.
  export_cols <- setdiff(names(df), c("x_centre", "x_pos"))
  export_data <- df[, export_cols, drop = FALSE]
  rownames(export_data) <- NULL

  if (!is.null(excel_filename)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      warning("openxlsx is not installed; skipping Excel data save.")
    } else {
      wb <- openxlsx::createWorkbook()
      openxlsx::modifyBaseFont(wb, fontSize = 12, fontColour = "black",
                               fontName = "Arial")
      openxlsx::addWorksheet(wb, "data")
      openxlsx::writeData(wb, "data", export_data)
      if (ncol(export_data) >= 1) {
        openxlsx::setColWidths(wb, "data",
                               cols = seq_len(ncol(export_data)),
                               widths = "auto")
        openxlsx::freezePane(wb, "data", firstRow = TRUE)
        openxlsx::addFilter(wb, "data", rows = 1,
                            cols = seq_len(ncol(export_data)))
      }
      openxlsx::saveWorkbook(wb, excel_filename, overwrite = TRUE)
      message("Saved data to: ", excel_filename)
    }
  }
  if (!is.null(output_file) && !is.null(title)) {
    message("Saved figure to: ", output_file)
  }

  # ── 14. Return base_plot object ──────────────────────────────────────────
  invisible(structure(
    list(draw = draw_fn, width = fig_w, height = fig_h, data = export_data),
    class = "base_plot"
  ))
}

# Allow access to CATALYST internal helpers (.check_sce, .check_k)
