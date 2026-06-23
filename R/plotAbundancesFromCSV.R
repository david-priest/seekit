# plotAbundancesFromCSV.R
#
# Thin CSV entry point used by the Shiny tuner.
# Requires that plotAbundancesBase.R has already been sourced.

plotAbundancesFromCSV <- function(
  points_csv,
  output_file      = NULL,
  svg_grouped      = TRUE,
  dodge_w          = 10.0,
  group_gap        = 2.0,
  jitter_w         = 0.5,
  n_cols           = 2L,
  per_panel_width  = 11.0,
  per_panel_height = 3.5,
  base_font_size   = 9.0,
  point_size       = 35.0,
  point_alpha      = 0.70,
  line_alpha       = 0.35,
  line_width       = 0.55,
  split_val_filter = NULL,
  group_order      = NULL,
  clusters_order   = NULL,
  patient_shapes   = NULL,   # named integer vector (patient -> pch); NULL = auto
  cmv_gap          = 0       # extra space inserted between CMV- and CMV+ conditions
) {
  if (!exists(".plotFromDF", mode = "function")) {
    stop(".plotFromDF not found. Source plotAbundancesBase.R before plotAbundancesFromCSV.R")
  }

  pts_all <- read.csv(
    as.character(points_csv),
    stringsAsFactors = FALSE,
    colClasses = c(cluster_id = "character", patient_id = "character")
  )

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
    clusters_order   = clusters_order,
    patient_shapes   = patient_shapes,
    cmv_gap          = cmv_gap
  )
}
