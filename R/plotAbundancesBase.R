# plotAbundancesBase.R
#
# Base-R strip-plot for cluster abundances — drop-in replacement for
# plotAbundances22 that uses base R graphics instead of ggplot2.
#
# Identical data-extraction logic to plotAbundances22 (CATALYST frequencies,
# same metadata join, same condition & eTreg handling, same split_col support).
# Identical layout maths to plot_abundances_csv.py.
# Requires CATALYST, SingleCellExperiment, dplyr.
#
# ─────────────────────────────────────────────────────────────────────────────
# TYPICAL CALL
# ─────────────────────────────────────────────────────────────────────────────
#
#   # Main-text figure only
#   plotAbundancesBase(
#     x            = sce_plot,
#     k            = "meta20",
#     group_by     = "T cell",
#     split_col    = "B cell",
#     meta_df      = meta_df,
#     etreg_donors = c("LP1", "LP9", "LP10"),
#     dodge_w      = 10,
#     group_gap    = 2,
#     jitter_w     = 0.5
#   )
#
#   # Main-text + supplementary replicates overlaid on the same figure.
#   # Replicate points share the same fill colours as main-text points but are
#   # drawn with a dark contrasting border so they are immediately recognisable.
#   plotAbundancesBase(
#     x            = sce_plot,          # plot1_mainplotmem samples
#     rep_sce      = sce_plot4,         # plot4_repplotmem samples
#     k            = "meta20",
#     group_by     = "T cell",
#     split_col    = "B cell",
#     meta_df      = meta_df,
#     etreg_donors = c("LP1", "LP9", "LP10"),
#     output_file  = "abundances_with_reps.pdf",
#     dodge_w      = 10,
#     group_gap    = 2,
#     jitter_w     = 0.5
#   )
#
# ─────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(SingleCellExperiment)
# svglite is referenced fully-qualified (svglite::svgstring / svglite::svglite)
# and is an optional (Suggests) dependency, so it is not attached here.

# ── Visual constants ──────────────────────────────────────────────────────────

COND_COLORS_BASE <- c(
  "DMSO CMV-" = "#999999",
  "pp65 CMV-" = "#4477BB",
  "DMSO CMV+" = "#999999",
  "pp65 CMV+" = "#CC4444",
  "eTreg"     = "#EBCC2A",
  "ctfh"      = "#44AA88"   # add further conditions here as needed
)

COND_ORDER_BASE <- c("DMSO CMV-", "pp65 CMV-", "DMSO CMV+", "pp65 CMV+", "eTreg", "ctfh")

SEGMENT_PAIRS_BASE <- list(
  c("DMSO CMV-", "pp65 CMV-"),
  c("DMSO CMV+", "pp65 CMV+"),
  c("pp65 CMV+", "eTreg")
)

# pch 21-25: filled shapes with separate border/fill (col = border, bg = fill)
SHAPE_CYCLE_BASE <- c(21L, 22L, 23L, 24L, 25L, 21L, 22L, 23L)

# Border colour used for replicate points (fill is the same as main-text points)
REP_BORDER_COLOR <- "#333333"
REP_BORDER_LWD   <- 1.5     # slightly heavier than the main-text border


# ══════════════════════════════════════════════════════════════════════════════
# Internal layout helpers  (identical maths to plot_abundances_csv.py)
# ══════════════════════════════════════════════════════════════════════════════

.computeLayout <- function(x_labels, group_gap, dodge_w, cond_levels,
                            strip_margin = 0.1, cmv_gap = 0) {
  n_cond  <- length(cond_levels)
  spacing <- dodge_w / n_cond

  # Evenly-spaced positions, then open an extra gap after the 2nd condition
  # (CMV− / CMV+ boundary).  When cmv_gap = 0 the result is identical to the
  # original formula.
  raw_pos <- (seq_len(n_cond) - 1L) * spacing
  if (cmv_gap > 0 && n_cond > 2L)
    raw_pos[-(1:2)] <- raw_pos[-(1:2)] + cmv_gap
  raw_offsets <- raw_pos - mean(raw_pos)
  pd_offsets  <- setNames(raw_offsets, cond_levels)

  # Strip width: at least as wide as dodge_w/2 (original behaviour), expanded
  # only if the gap pushes offsets beyond that.
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

.addXPos <- function(df, centers, pd_offsets, jitter_w = 0, seed = 42L) {
  df$x_centre <- unname(centers[as.character(df$x_label)])
  df$x_pos    <- df$x_centre + unname(pd_offsets[as.character(df$condition)])
  if (jitter_w > 0) {
    set.seed(seed)
    df$x_pos <- df$x_pos + runif(nrow(df), -jitter_w / 2, jitter_w / 2)
  }
  df
}

.buildSegments <- function(pts) {
  join_keys <- c("cluster_id", "patient_id", "x_label")
  if ("exp"       %in% names(pts)) join_keys <- c(join_keys, "exp")
  if ("split_val" %in% names(pts)) join_keys <- c(join_keys, "split_val")
  # Keep main and replicate segments separate so lines never cross groups
  if ("is_rep"    %in% names(pts)) join_keys <- c(join_keys, "is_rep")

  rows <- lapply(SEGMENT_PAIRS_BASE, function(pair) {
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
# SVG composition helpers
#
# Each panel (cluster plot + legend) is rendered into its own in-memory SVG
# via svglite::svgstring, then stitched together into a single output SVG
# where every panel lives inside its own <g id="..."> group.  This gives
# Illustrator / Inkscape users cleanly separable panels.
# ══════════════════════════════════════════════════════════════════════════════

.renderPanelSVG <- function(draw_fn, w_in, h_in, pointsize = 12) {
  s <- svglite::svgstring(width = w_in, height = h_in, pointsize = pointsize,
                          standalone = FALSE)
  done <- FALSE
  on.exit(if (!done) try(dev.off(), silent = TRUE))
  draw_fn()
  dev.off()
  done <- TRUE
  as.character(s())
}

.extractSVGBody <- function(svg_str) {
  svg_start <- regexpr("<svg", svg_str, fixed = TRUE)
  if (svg_start < 0) return(list(inner = "", width = NA, height = NA))
  tail <- substring(svg_str, svg_start)
  tag_end <- regexpr(">", tail, fixed = TRUE)
  open_tag <- substr(tail, 1L, tag_end)
  rest <- substring(tail, tag_end + 1L)
  close_pos <- regexpr("</svg>", rest, fixed = TRUE)
  inner <- if (close_pos > 0) substr(rest, 1L, close_pos - 1L) else rest

  get_num <- function(attr) {
    m <- regmatches(open_tag,
                    regexpr(paste0(attr, '="[0-9.]+'), open_tag))
    if (length(m) == 0L) return(NA_real_)
    as.numeric(sub(paste0(attr, '="'), "", m, fixed = TRUE))
  }
  list(inner = inner, width = get_num("width"), height = get_num("height"))
}

.namespaceSVGIds <- function(inner, prefix) {
  # Extract id values from both double- and single-quoted attributes.
  m_dq <- regmatches(inner, gregexpr('id="[^"]+"', inner, perl = TRUE))[[1L]]
  m_sq <- regmatches(inner, gregexpr("id='[^']+'", inner, perl = TRUE))[[1L]]

  ids <- character()
  if (length(m_dq) > 0L && !identical(m_dq, ""))
    ids <- c(ids, sub('^id="([^"]+)"$', '\\1', m_dq, perl = TRUE))
  if (length(m_sq) > 0L && !identical(m_sq, ""))
    ids <- c(ids, sub("^id='([^']+)'$", "\\1", m_sq, perl = TRUE))

  ids <- unique(ids[nzchar(ids)])
  if (length(ids) == 0L) return(inner)

  prefix <- .safeId(prefix)
  for (old_id in ids) {
    new_id <- paste0(prefix, "__", old_id)

    inner <- gsub(paste0('id="', old_id, '"'),
                  paste0('id="', new_id, '"'), inner, fixed = TRUE)
    inner <- gsub(paste0("id='", old_id, "'"),
                  paste0("id='", new_id, "'"), inner, fixed = TRUE)

    # Update common SVG references to the renamed ids.
    inner <- gsub(paste0("url(#", old_id, ")"),
                  paste0("url(#", new_id, ")"), inner, fixed = TRUE)
    inner <- gsub(paste0('href="#', old_id, '"'),
                  paste0('href="#', new_id, '"'), inner, fixed = TRUE)
    inner <- gsub(paste0("href='#", old_id, "'"),
                  paste0("href='#", new_id, "'"), inner, fixed = TRUE)
    inner <- gsub(paste0('xlink:href="#', old_id, '"'),
                  paste0('xlink:href="#', new_id, '"'), inner, fixed = TRUE)
    inner <- gsub(paste0("xlink:href='#", old_id, "'"),
                  paste0("xlink:href='#", new_id, "'"), inner, fixed = TRUE)
  }

  inner
}

.removeCanvasClip <- function(inner) {
  # svglite emits a full-device clipPath (id commonly containing "cpMC4...")
  # that can clip rotated axis labels after panel stitching. Remove only that
  # outer clip while keeping plot-region clipPaths intact.
  inner <- gsub(
    "<defs>\\s*<clipPath id='[^']*cpMC4[^']*'>\\s*<rect x='0\\.00' y='0\\.00' width='[0-9.]+' height='[0-9.]+' ?/>\\s*</clipPath>\\s*</defs>",
    "",
    inner,
    perl = TRUE
  )
  inner <- gsub(
    '<defs>\\s*<clipPath id="[^"]*cpMC4[^"]*">\\s*<rect x="0\\.00" y="0\\.00" width="[0-9.]+" height="[0-9.]+" ?/>\\s*</clipPath>\\s*</defs>',
    "",
    inner,
    perl = TRUE
  )

  inner <- gsub("clip-path='url\\(#[^']*cpMC4[^']*\\)'", "", inner, perl = TRUE)
  inner <- gsub('clip-path="url\\(#[^"]*cpMC4[^"]*\\)"', "", inner, perl = TRUE)
  inner
}

.removePanelWhiteBackground <- function(inner) {
  # Remove svglite full-panel white background rectangles so stitched panels
  # remain transparent unless the caller adds an explicit background.
  inner <- gsub(
    "<rect\\s+width='100%'\\s+height='100%'\\s+style='[^']*fill:\\s*#FFFFFF;?[^']*'\\s*/>",
    "",
    inner,
    perl = TRUE
  )
  inner <- gsub(
    '<rect\\s+width="100%"\\s+height="100%"\\s+style="[^"]*fill:\\s*#FFFFFF;?[^"]*"\\s*/>',
    "",
    inner,
    perl = TRUE
  )
  inner
}

.xmlEscape <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;",  s, fixed = TRUE)
  s <- gsub(">", "&gt;",  s, fixed = TRUE)
  s
}

.safeId <- function(s) {
  s <- gsub("[^A-Za-z0-9_-]+", "_", as.character(s))
  if (!nzchar(s)) "panel" else s
}

.composeGroupedSVG <- function(sections, n_cols, pw_pt, ph_pt, output_file,
                               panel_gap_x_pt = 20,
                               panel_gap_y_pt = 16,
                               section_gap_pt = 12,
                               panel_pad_left_pt = 24,
                               panel_pad_right_pt = 12,
                               panel_pad_top_pt = 16,
                               panel_pad_bottom_pt = 30) {
  # sections: list of list(title = <char|NULL>, panels = list(list(id, inner)))
  # pw_pt / ph_pt : panel width / height in SVG user units (pt, i.e. inches*72)
  title_h <- 22
  holder_w <- pw_pt + panel_pad_left_pt + panel_pad_right_pt
  holder_h <- ph_pt + panel_pad_top_pt + panel_pad_bottom_pt
  total_w <- holder_w * n_cols + panel_gap_x_pt * max(0L, n_cols - 1L)
  y_cursor <- 0
  body <- character()

  for (s_idx in seq_along(sections)) {
    sec <- sections[[s_idx]]
    if (!is.null(sec$title) && nzchar(sec$title)) {
      body <- c(body, sprintf(
        '<g id="%s"><text x="%g" y="%g" text-anchor="middle" font-family="Arial" font-weight="bold" font-size="14">%s</text></g>',
        .safeId(paste0("title_", sec$title)),
        total_w / 2, y_cursor + title_h * 0.75,
        .xmlEscape(sec$title)))
      y_cursor <- y_cursor + title_h
    }

    n <- length(sec$panels)
    n_rows <- ceiling(n / n_cols)
    for (i in seq_len(n)) {
      col <- (i - 1L) %% n_cols
      row <- (i - 1L) %/% n_cols
      x   <- col * (holder_w + panel_gap_x_pt)
      y   <- y_cursor + row * (holder_h + panel_gap_y_pt)
      p   <- sec$panels[[i]]
      body <- c(body,
        sprintf('<g id="%s" transform="translate(%g, %g)">', .safeId(p$id), x, y),
        sprintf('<svg class="svglite" x="0" y="0" width="%g" height="%g" viewBox="-%g -%g %g %g" overflow="hidden">',
                holder_w, holder_h,
                panel_pad_left_pt, panel_pad_top_pt,
                holder_w, holder_h),
        p$inner,
        '</svg>',
        '</g>')
    }
    y_cursor <- y_cursor + n_rows * holder_h + max(0L, n_rows - 1L) * panel_gap_y_pt
    if (s_idx < length(sections))
      y_cursor <- y_cursor + section_gap_pt
  }

  total_h <- y_cursor
  header <- c(
    '<?xml version="1.0" encoding="UTF-8"?>',
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="%gpt" height="%gpt" viewBox="0 0 %g %g">',
            total_w, total_h, total_w, total_h))
  writeLines(c(header, body, '</svg>'), output_file)
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal drawing helpers
# ══════════════════════════════════════════════════════════════════════════════

# Draw a margin label wrapped in a TIGHT clip rectangle around its bounding
# box. The default base-R behaviour for `title()` / `mtext()` is to emit
# margin text inside a clip path the size of the figure region (the mfrow
# cell). When the PDF opens in Illustrator that translates to a panel-sized
# clipping mask around every label — annoying to edit. This helper issues
# an explicit `clip()` rectangle just before drawing, so the per-label clip
# in the PDF is the text's own bounding box plus a small pad.
#
# After drawing, clip() is restored to the plot region so subsequent draw
# calls (axes, points, etc.) aren't constrained.
#
# Args:
#   label   character to draw
#   at      x position in user coords (default: midpoint of x usr range)
#   side    margin side: 1 = bottom, 3 = top (matches mtext)
#   line    margin-line offset from the plot edge (matches mtext: positive
#           = away from plot, negative = inside plot from that edge)
#   cex     char expansion
#   font    1 = plain, 2 = bold, etc.
.drawClippedMarginText <- function(label, at = NULL, side = 3,
                                    line = 1.4, cex = 1, font = 1) {
  usr <- par("usr")
  pin <- par("pin")          # plot region size in inches  (w, h)
  cin <- par("cin")          # default char size in inches (w, h)
  if (is.null(at)) at <- mean(usr[1:2])

  # 1 mtext "line" = cin[2] inches vertically. Convert to user-y units.
  lh_user_y <- cin[2] * (diff(usr[3:4]) / pin[2])
  y_user <- switch(as.character(side),
                   "1" = usr[3] - line * lh_user_y,
                   "3" = usr[4] + line * lh_user_y,
                   stop(".drawClippedMarginText: side must be 1 (bottom) or 3 (top)"))

  tw <- strwidth(label, cex = cex, font = font)
  th <- strheight(label, cex = cex, font = font)
  pad_x <- tw * 0.06     # ~6% horizontal pad each side
  pad_y <- th * 0.20     # ~20% vertical pad (text-baseline asymmetry)

  # Allow drawing into the margin (xpd = TRUE) then narrow the clip to the
  # text's own bounding box. `clip()` is the authoritative clip rect; xpd
  # just controls whether margin-drawing is permitted at all.
  old_xpd <- par("xpd")
  par(xpd = TRUE)
  on.exit(par(xpd = old_xpd), add = TRUE)

  clip(at - tw/2 - pad_x, at + tw/2 + pad_x,
       y_user - th/2 - pad_y, y_user + th/2 + pad_y)
  text(at, y_user, labels = label, cex = cex, font = font)

  # Restore clip to the plot region for subsequent draw calls.
  clip(usr[1], usr[2], usr[3], usr[4])
}


.drawFacet <- function(c_pts, c_segs, layout_df, cluster,
                        shape_map, x_lim,
                        point_size, point_alpha, line_alpha, line_width,
                        font_cex, pd_offsets, cond_levels,
                        title_scale = 1.0, group_label_scale = 0.85,
                        y_axis_title = NULL) {

  y_max <- if (nrow(c_pts) > 0L) max(c_pts$Freq, na.rm = TRUE) else 1.0
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1.0
  y_pad <- y_max * 0.08
  y_lim <- c(-y_pad, y_max * 1.18 + 1e-4)

  graphics::plot.new()
  graphics::plot.window(xlim = x_lim, ylim = y_lim)

  # Faint horizontal gridlines at the major y-ticks (drawn FIRST so that
  # segments + points render on top). Helps the eye read absolute %
  # values when n is small and the data is sparse.
  abline(h = axTicks(2), col = "#eeeeee", lwd = 0.3)

  # Dashed vertical separators between T cell groups
  if (nrow(layout_df) > 1L) {
    for (i in seq_len(nrow(layout_df) - 1L)) {
      x_sep <- (layout_df$strip_xmax[i] + layout_df$strip_xmin[i + 1L]) / 2
      abline(v = x_sep, col = "#bbbbbb", lty = 2L, lwd = 0.8)
    }
  }

  # Connecting segments — replicate segments drawn dashed so they are visually
  # distinct even before looking at the point borders.
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
               col = adjustcolor(REP_BORDER_COLOR, alpha.f = line_alpha),
               lwd = line_width * 2.0, lty = 2L)
  }

  # Points
  # pch 21-25: filled shapes — col = border colour, bg = fill colour.
  # Border = same colour as fill but FULL opacity; bg = same colour with
  # `point_alpha` applied. Gives every point a crisp coloured edge against
  # an alphad interior, which keeps overlapping points distinguishable
  # without the previous "white halo" look. Replicate points keep their
  # special REP_BORDER_COLOR outline so they remain visually distinct.
  if (nrow(c_pts) > 0L) {
    pt_cex <- sqrt(point_size / 20)
    for (i in seq_len(nrow(c_pts))) {
      r      <- c_pts[i, ]
      fill   <- COND_COLORS_BASE[as.character(r$condition)]
      if (is.na(fill)) fill <- "#888888"
      pch    <- shape_map[[as.character(r$patient_id)]]
      if (is.null(pch) || is.na(pch)) pch <- 21L
      is_rep      <- isTRUE(r$is_rep)
      is_fillable <- pch %in% 21:25
      border <- if (is_fillable) {
        if (is_rep) REP_BORDER_COLOR else fill  # was "white"
      } else {
        fill
      }
      lwd <- if (is_fillable) {
        if (is_rep) REP_BORDER_LWD else 0.7   # slightly thicker for the new
                                              # same-colour border to read
      } else 1.5
      points(r$x_pos, r$Freq,
             pch = pch, cex = pt_cex,
             col = border,  # full-opacity outline (was alpha'd)
             bg  = adjustcolor(if (is_fillable) fill else NA, alpha.f = point_alpha),
             lwd = lwd)
    }
  }

  # Multipliers chosen to be ≥0.85 (was 0.65) so that at typical paper-figure
  # downscale (~6× from rendering width to a 2-column figure) the smallest
  # text still lands above 6 pt. Bump base_font_size in the chunk if you
  # need uniformly larger labels. y-axis title drawn ONLY when
  # y_axis_title is non-NULL (default off — one title per panel is usually
  # redundant once the y-tick numbers are visible).
  axis(2, cex.axis = font_cex * 0.9, las = 1, lwd = 0.7)
  if (!is.null(y_axis_title) && nzchar(y_axis_title)) {
    mtext(y_axis_title, side = 2, line = 2.8, cex = font_cex * 0.9)
  }

  # Build condition x positions: one per (T cell group × condition)
  cond_x   <- unlist(lapply(layout_df$x_centre,
                            function(xc) xc + unname(pd_offsets[cond_levels])))
  cond_lab <- rep(cond_levels, nrow(layout_df))

  # Small ticks at each condition position only
  axis(1, at = cond_x, labels = FALSE, tck = -0.015, lwd = 0.5)

  # Bottom-margin label positioning is line-based (constant inches below
  # the axis) so it stays correct when per_panel_height / margins change.
  # See plotAbundancesSimple.R for the rationale.
  usr <- par("usr")
  pin <- par("pin")
  cin <- par("cin")
  lh_user_y <- cin[2] * (diff(usr[3:4]) / pin[2])

  # Condition labels at bottom (45°). 0.5 lines below the axis; the label
  # then projects down-left ~5-6 lines depending on length.
  # adj = c(1, 0.5) keeps the rotated label right-anchored at the tick.
  text(cond_x, usr[3] - 0.5 * lh_user_y,
       labels = cond_lab,
       srt = 45, adj = c(1, 0.5), xpd = NA, cex = font_cex * 0.9)

  box(bty = "l", lwd = 0.7)

  # Panel title. Drawn via .drawClippedMarginText so it sits in a tight
  # per-label PDF clip rect (Illustrator shows a clip mask the size of
  # the title text, not the whole panel). line = 1.4 so it sits clear of
  # the T-cell group sub-labels at line = -0.2.
  .drawClippedMarginText(cluster, side = 3, line = 1.4,
                         cex = font_cex * title_scale, font = 2)

  # T cell group labels just below the panel title, inside the top margin
  for (i in seq_len(nrow(layout_df))) {
    mtext(layout_df$x_label[i], side = 3, at = layout_df$x_centre[i],
          line = -0.2, cex = font_cex * group_label_scale,
          font = 2L, xpd = NA)
  }
}


.drawLegend <- function(present_conds, present_patients, shape_map, font_cex,
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

  # Condition entries: white border + colour fill (main-text style)
  leg_col   <- c(rep("white", nc), rep("gray40", np))
  leg_pt.bg <- c(unname(COND_COLORS_BASE[present_conds]), rep(NA_character_, np))
  leg_lwd   <- c(rep(0.5, nc), rep(NA, np))
  leg_lab   <- c(present_conds, sorted_pats)

  if (has_rep) {
    # Two extra entries illustrating the main vs replicate border convention
    leg_lab   <- c(leg_lab,   "— Main text data",    "-- Replicate data")
    leg_pch   <- c(leg_pch,   22L,                    22L)
    leg_col   <- c(leg_col,   "white",                REP_BORDER_COLOR)
    leg_pt.bg <- c(leg_pt.bg, "gray60",               "gray60")
    leg_lwd   <- c(leg_lwd,   0.5,                    REP_BORDER_LWD)
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
# Internal: SCE → data frame (shared by main and rep pipelines)
# ══════════════════════════════════════════════════════════════════════════════
#
# Returns a data frame with columns:
#   cluster_id, sample_id, Freq, patient_id, x_label, split_val, condition
#   + any extra columns from metadata (exp, etc.)

.processSCE <- function(
  x,
  k,
  by,
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
  # 1. Frequencies
  c_ids <- .wl_cluster_ids(x, k)
  s_ids <- .wl_sample_ids(x)
  ns    <- table(cluster_id = c_ids, sample_id = s_ids)
  fq    <- prop.table(ns, 2) * 100
  df    <- as.data.frame(fq, stringsAsFactors = FALSE) %>%
    dplyr::filter(!is.nan(Freq))

  # 2. Metadata join
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

  # 3. Condition variable
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

    # Rows that did not match any of the above (e.g. ctfh or other conditions
    # present only in the replicate SCE).  Pass them through using the raw
    # stim value as the condition label so they appear on the plot.
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

  # 4. x_label column
  if (!group_by %in% names(df))
    stop("group_by column '", group_by, "' not found in data after metadata join.")
  df$x_label <- as.character(df[[group_by]])

  # 5. split_val column
  if (!is.null(split_col) && split_col %in% names(df)) {
    df$split_val <- as.character(df[[split_col]])
  } else {
    df$split_val <- "all"
  }

  # 6. Normalise patient_id / cluster_id
  if (!is.null(shape_by) && shape_by %in% names(df) && shape_by != "patient_id")
    df$patient_id <- as.character(df[[shape_by]])
  df$patient_id <- as.character(df$patient_id)
  df$cluster_id <- as.character(df$cluster_id)

  df
}


# ══════════════════════════════════════════════════════════════════════════════
# Internal: shared plotting engine (called by both SCE and CSV entry points)
# ══════════════════════════════════════════════════════════════════════════════

.plotFromDF <- function(
  pts_all,            # data frame: cluster_id, patient_id, x_label,
                      #   condition, Freq, split_val  [+ optional exp, is_rep]
  output_file       = NULL,
  svg_grouped       = TRUE,
  dodge_w           = 10.0,
  group_gap         = 2.0,
  jitter_w          = 0.5,
  n_cols            = 2L,
  per_panel_width   = 11.0,
  per_panel_height  = 3.5,
  base_font_size    = 9.0,
  point_size        = 35.0,
  point_alpha       = 0.70,
  line_alpha        = 0.35,
  line_width        = 0.55,
  split_val_filter  = NULL,
  group_order       = NULL,
  clusters_order    = NULL,
  patient_shapes    = NULL,  # named integer vector (patient → pch); NULL = auto
  cmv_gap           = 0,     # extra space inserted between CMV− and CMV+ conditions
  title_scale       = 1.0,   # cex multiplier for panel titles (relative to base_font_size)
  group_label_scale = 0.85,  # cex multiplier for T cell group labels
  y_axis_title      = NULL   # per-panel y-axis title (NULL = off, default)
) {
  pts_all$cluster_id <- as.character(pts_all$cluster_id)
  pts_all$patient_id <- as.character(pts_all$patient_id)
  if (!"is_rep" %in% names(pts_all)) pts_all$is_rep <- FALSE

  has_rep <- any(pts_all$is_rep, na.rm = TRUE)

  cond_levels    <- intersect(COND_ORDER_BASE, unique(pts_all$condition))
  all_x_labels   <- sort(unique(as.character(pts_all$x_label)))
  ordered_labels <- if (!is.null(group_order)) {
    c(intersect(group_order, all_x_labels), setdiff(all_x_labels, group_order))
  } else all_x_labels

  lay        <- .computeLayout(ordered_labels, group_gap, dodge_w, cond_levels,
                               cmv_gap = cmv_gap)
  centers    <- lay$centers
  pd_offsets <- lay$pd_offsets
  layout_df  <- lay$layout_df

  x_margin <- 0.5
  x_lim    <- c(min(layout_df$strip_xmin) - x_margin,
                max(layout_df$strip_xmax) + x_margin)

  pts_all <- .addXPos(pts_all, centers, pd_offsets, jitter_w = jitter_w)

  all_patients <- sort(unique(pts_all$patient_id))
  if (is.null(patient_shapes)) {
    shape_map <- setNames(
      as.list(SHAPE_CYCLE_BASE[
        (seq_along(all_patients) - 1L) %% length(SHAPE_CYCLE_BASE) + 1L
      ]),
      all_patients
    )
  } else {
    shape_map <- as.list(patient_shapes)
    missing_p <- setdiff(all_patients, names(shape_map))
    for (i in seq_along(missing_p))
      shape_map[[missing_p[i]]] <-
        SHAPE_CYCLE_BASE[(i - 1L) %% length(SHAPE_CYCLE_BASE) + 1L]
  }

  all_clusters     <- sort(unique(pts_all$cluster_id))
  ordered_clusters <- if (!is.null(clusters_order)) {
    c(intersect(clusters_order, all_clusters), setdiff(all_clusters, clusters_order))
  } else all_clusters

  split_vals <- sort(unique(as.character(pts_all$split_val)))
  if (!is.null(split_val_filter))
    split_vals <- intersect(split_vals, split_val_filter)

  font_cex <- base_font_size / 12

  # Calculate figure dimensions from per-panel sizes
  n_cl  <- length(ordered_clusters)
  n_rows_fig <- ceiling((n_cl + 1L) / n_cols)
  fig_w <- per_panel_width  * n_cols + 1.5
  fig_h <- per_panel_height * n_rows_fig

  drawSplitToCurrentDevice <- function(sv_pts, sv, sv_segs, present_clusters) {
    n_cells <- length(present_clusters) + 1L
    n_rows  <- ceiling(n_cells / n_cols)

    # Margins (in mtext-style "lines" — 1 line ≈ cin[2] inches ≈ 14.4 pt at
    # the device's default 12pt size). Total vertical gap between adjacent
    # panels in mfrow = mar[1] of upper + mar[3] of lower.
    #
    #   bottom: 8.0   — 45° rotated condition labels at cex 0.9 × font_cex
    #                   project ~5-6 line-heights below the plot region;
    #                   8 lines gives a buffer so they don't crash into the
    #                   panel below. If you increase base_font_size beyond
    #                   ~22, bump this proportionally.
    #   left  : 4.5   — room for the leftmost rotated label's tail
    #                   ("DMSO CMV−") and the y-axis tick numbers
    #   top   : 5.5   — fits the panel title (line = 1.4, cex up to ~1.5)
    #                   AND the T-cell group sub-labels at line = -0.2,
    #                   with breathing room above the previous panel's
    #                   bottom labels
    #   right : 0.5   — padding so the rightmost "eTreg" label doesn't
    #                   clip at the panel edge
    #
    # If panels still feel cramped, the most effective lever is the chunk's
    # `per_panel_height` argument — bump from 4 to 5–6 inches and the plot
    # region grows back into the increased margin space.
    old_par <- par(
      mfrow = c(n_rows, n_cols),
      mar   = c(8.0, 4.5, 5.5, 0.5),
      oma   = c(0, 0, if (sv != "all") 1.8 else 0, 0)
    )
    on.exit(par(old_par), add = TRUE)

    for (cluster in present_clusters) {
      c_pts  <- sv_pts[sv_pts$cluster_id == cluster, , drop = FALSE]
      c_segs <- if (nrow(sv_segs) > 0L)
                  sv_segs[sv_segs$cluster_id == cluster, , drop = FALSE]
                else data.frame()
      .drawFacet(c_pts, c_segs, layout_df, cluster,
                 shape_map, x_lim,
                 point_size, point_alpha, line_alpha, line_width, font_cex,
                 pd_offsets, cond_levels,
                 title_scale = title_scale, group_label_scale = group_label_scale,
                 y_axis_title = y_axis_title)
    }

    present_conds    <- intersect(COND_ORDER_BASE, unique(sv_pts$condition))
    present_patients <- sort(unique(sv_pts$patient_id))
    .drawLegend(present_conds, present_patients, shape_map, font_cex,
                has_rep = has_rep)

    n_empty <- n_rows * n_cols - n_cells
    if (n_empty > 0L) {
      for (i in seq_len(n_empty)) {
        par(mar = c(0, 0, 0, 0))
        plot.new()
      }
    }

    if (sv != "all")
      mtext(sv, outer = TRUE, cex = font_cex * 1.3, font = 2, line = 1)
  }

  # ── Path A: no output_file → draw to current device (screen / existing dev)
  if (is.null(output_file)) {
    for (sv in split_vals) {
      sv_pts <- pts_all[as.character(pts_all$split_val) == sv, , drop = FALSE]
      if (nrow(sv_pts) == 0L) next
      present_clusters <- intersect(ordered_clusters, unique(sv_pts$cluster_id))
      if (length(present_clusters) == 0L) next

      sv_segs  <- .buildSegments(sv_pts)
      drawSplitToCurrentDevice(sv_pts, sv, sv_segs, present_clusters)
    }
    return(invisible(recordPlot()))
  }

  output_ext <- tolower(tools::file_ext(output_file))

  # ── Path B: standard graphics device layout (including ungrouped SVG) ─────
  if (output_ext != "svg" || !isTRUE(svg_grouped)) {
    if (output_ext == "pdf") {
      grDevices::pdf(output_file, width = fig_w, height = fig_h, onefile = TRUE)
    } else if (output_ext == "svg") {
      svglite::svglite(output_file, width = fig_w, height = fig_h, pointsize = 12)
    } else if (output_ext == "png") {
      grDevices::png(output_file, width = fig_w, height = fig_h,
                     units = "in", res = 300)
    } else if (output_ext %in% c("jpg", "jpeg")) {
      grDevices::jpeg(output_file, width = fig_w, height = fig_h,
                      units = "in", res = 300)
    } else if (output_ext %in% c("tif", "tiff")) {
      grDevices::tiff(output_file, width = fig_w, height = fig_h,
                      units = "in", res = 300)
    } else if (output_ext == "bmp") {
      grDevices::bmp(output_file, width = fig_w, height = fig_h,
                     units = "in", res = 300)
    } else {
       stop("Unsupported output format: '", output_ext,
         "'. Use .svg/.pdf or a standard raster extension.")
    }
    on.exit(dev.off(), add = TRUE)

    for (sv in split_vals) {
      sv_pts <- pts_all[as.character(pts_all$split_val) == sv, , drop = FALSE]
      if (nrow(sv_pts) == 0L) next
      present_clusters <- intersect(ordered_clusters, unique(sv_pts$cluster_id))
      if (length(present_clusters) == 0L) next

      sv_segs <- .buildSegments(sv_pts)
      drawSplitToCurrentDevice(sv_pts, sv, sv_segs, present_clusters)
    }
    return(invisible(NULL))
  }

  # ── Path C: grouped-SVG output (one <g id="..."> per panel) ───────────────
  render_one <- function(draw_fn) {
    .renderPanelSVG(draw_fn,
                    w_in = per_panel_width,
                    h_in = per_panel_height,
                    pointsize = 12)
  }

  sections <- list()

  for (sv in split_vals) {
    sv_pts <- pts_all[as.character(pts_all$split_val) == sv, , drop = FALSE]
    if (nrow(sv_pts) == 0L) next
    present_clusters <- intersect(ordered_clusters, unique(sv_pts$cluster_id))
    if (length(present_clusters) == 0L) next

    sv_segs <- .buildSegments(sv_pts)
    panels  <- list()

    for (cluster in present_clusters) {
      c_pts  <- sv_pts[sv_pts$cluster_id == cluster, , drop = FALSE]
      c_segs <- if (nrow(sv_segs) > 0L)
                  sv_segs[sv_segs$cluster_id == cluster, , drop = FALSE]
                else data.frame()

      svg_str <- render_one(function() {
        par(mar = c(4.5, 3.0, 2.2, 0.3))
        .drawFacet(c_pts, c_segs, layout_df, cluster,
                   shape_map, x_lim,
                   point_size, point_alpha, line_alpha, line_width, font_cex,
                   pd_offsets, cond_levels,
                   title_scale = title_scale,
                   group_label_scale = group_label_scale,
                   y_axis_title = y_axis_title)
      })
      body <- .extractSVGBody(svg_str)
      panel_id <- paste0(sv, "__", cluster)
      panel_inner <- .namespaceSVGIds(body$inner, paste0("panel_", panel_id))
      panel_inner <- .removeCanvasClip(panel_inner)
      panel_inner <- .removePanelWhiteBackground(panel_inner)
      panels[[length(panels) + 1L]] <-
        list(id = panel_id, inner = panel_inner)
    }

    present_conds    <- intersect(COND_ORDER_BASE, unique(sv_pts$condition))
    present_patients <- sort(unique(sv_pts$patient_id))
    legend_str <- render_one(function() {
      .drawLegend(present_conds, present_patients, shape_map, font_cex,
                  has_rep = has_rep)
    })
    body <- .extractSVGBody(legend_str)
    legend_id <- paste0(sv, "__legend")
    legend_inner <- .namespaceSVGIds(body$inner, paste0("panel_", legend_id))
    legend_inner <- .removeCanvasClip(legend_inner)
    legend_inner <- .removePanelWhiteBackground(legend_inner)
    panels[[length(panels) + 1L]] <-
      list(id = legend_id, inner = legend_inner)

    sections[[length(sections) + 1L]] <- list(
      title  = if (sv != "all") sv else NULL,
      panels = panels
    )
  }

  dpi <- 72   # svglite emits SVG user-units as inches * 72
  .composeGroupedSVG(
    sections,
    n_cols  = n_cols,
    pw_pt   = per_panel_width  * dpi,
    ph_pt   = per_panel_height * dpi,
    panel_gap_x_pt = 20,
    panel_gap_y_pt = 16,
    section_gap_pt = 12,
    output_file = output_file
  )

  invisible(NULL)
}


# ══════════════════════════════════════════════════════════════════════════════
# plotAbundancesBase  — main entry point, takes SCE directly
# ══════════════════════════════════════════════════════════════════════════════
#
# All parameters identical to plotAbundances22 plus:
#   rep_sce          optional SCE of technical-replicate samples (plot4).
#                    When supplied these are overlaid on every panel with the
#                    same condition fill colours but a dark contrasting border
#                    (REP_BORDER_COLOR) so they are visually distinct from the
#                    main-text points.  Connecting lines for replicate samples
#                    are drawn dashed.  The legend gains "Main text data" /
#                    "Replicate data" entries.
#   output_file      path for PDF/PNG output  (NULL = current device)
#   svg_grouped      when output_file is .svg, TRUE writes grouped per-panel
#                    SVG (best for selective editing), FALSE writes a plain
#                    ungrouped SVG directly from base graphics (often faster
#                    in Illustrator for quick edits).
#   jitter_w         sub-dodge jitter width    (default 0.5)
#   per_panel_width  panel width in inches      (default 11)
#   per_panel_height panel height in inches     (default 3.5)
#   group_gap        white space between strips (replaces group_spacing)

plotAbundancesBase <- function(
  x,
  rep_sce           = NULL,            # optional: replicate SCE (plot4 samples)
  k                 = "meta20",
  by                = c("sample_id", "cluster_id"),
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
  svg_grouped       = TRUE,
  per_panel_width   = 11.0,
  per_panel_height  = 3.5,
  base_font_size    = 9.0,
  point_size        = 35.0,
  point_alpha       = 0.70,
  split_val_filter  = NULL,
  patient_shapes    = NULL,            # named integer vector (patient → pch); NULL = auto
  cmv_gap           = 0,              # extra space inserted between CMV− and CMV+ conditions
  title_scale       = 1.0,            # cex multiplier for panel titles
  group_label_scale = 0.85,           # cex multiplier for T cell group labels
  # Per-panel y-axis title (e.g. "Proportion [%]"). NULL = no label, which
  # is the default — the y-tick numbers already convey scale and a label
  # per panel is usually visual noise. Pass a string to enable.
  y_axis_title      = NULL,
  # ---- f2-style auto-save (added 2026-05-25) -------------------------------
  # When `title` is non-NULL the function auto-computes an `output_file` in
  # the standard dated folder (here::here(paste0(Sys.Date(), "_", rmd_name)))
  # and saves the figure there, mirroring f2() / saveFig() conventions. The
  # default `format = "pdf"` produces an editable PDF (preferred over SVG
  # for Illustrator workflows). Pass `saveExcel = TRUE` to also write the
  # underlying long-format point dataframe to xlsx alongside.
  #
  # Explicit `output_file = ...` still wins (backwards-compatible). If both
  # `title` and `output_file` are set, a warning is emitted and output_file
  # takes precedence.
  title             = NULL,
  saveExcel         = FALSE,
  out_format        = c("pdf", "svg")
) {
  by <- match.arg(by)
  out_format <- match.arg(out_format)

  .wl_check_sce(x, TRUE)
  k <- .wl_check_k(x, k)

  use_etreg <- !is.null(etreg_donors)

  # ---- f2-style folder + filename construction ---------------------------
  # Compute output_file and excel_filename from `title` if provided. Uses
  # here::here() so the folder lands in the qmd's project root regardless
  # of the current working directory (matches the fix applied to f.R / f2.R
  # / fl2.R / fx.R / fxl.R — see ledger G13).
  excel_filename <- NULL
  if (!is.null(title)) {
    if (!is.null(output_file)) {
      warning("plotAbundancesBase: both `title` and `output_file` provided. ",
              "`output_file` takes precedence for the figure file; ",
              "`saveExcel` will use `title` for the data filename.")
    }
    if (!requireNamespace("rmdhelp", quietly = TRUE)) {
      rmd_name_noext <- "plotAbundancesBase"
    } else {
      rmd_file_path <- tryCatch(rmdhelp::get_this_rmd_file(),
                                error = function(e) NULL)
      rmd_name_noext <- if (!is.null(rmd_file_path))
        tools::file_path_sans_ext(basename(rmd_file_path)) else "plotAbundancesBase"
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

  # Shared argument list for .processSCE
  proc_args <- list(
    k            = k,
    by           = by,
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

  # ── Process main SCE ────────────────────────────────────────────────────────
  df_main        <- do.call(.processSCE, c(list(x = x), proc_args))
  df_main$is_rep <- FALSE

  # ── Process replicate SCE (optional) ────────────────────────────────────────
  if (!is.null(rep_sce)) {
    .wl_check_sce(rep_sce, TRUE)
    df_rep        <- do.call(.processSCE, c(list(x = rep_sce), proc_args))
    df_rep$is_rep <- TRUE
    pts_all       <- dplyr::bind_rows(df_main, df_rep)
  } else {
    pts_all <- df_main
  }

  # ── Draw ───────────────────────────────────────────────────────────────────
  .plotFromDF(
    pts_all          = pts_all,
    output_file      = output_file,
    svg_grouped      = svg_grouped,
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
    clusters_order    = clusters_order,
    patient_shapes    = patient_shapes,
    cmv_gap           = cmv_gap,
    title_scale       = title_scale,
    group_label_scale = group_label_scale,
    y_axis_title      = y_axis_title
  )

  # ── Save underlying data to Excel (if requested) ───────────────────────────
  # Mirrors f2()'s saveExcel = TRUE behaviour. Writes the long-format
  # per-point dataframe (one row per sample × cluster × stim × ...) to a
  # single "data" sheet. If a replicate SCE was provided, the `is_rep`
  # column distinguishes main vs replicate rows.
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
