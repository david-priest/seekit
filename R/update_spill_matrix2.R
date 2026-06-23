#' @title Manually override a spillover-matrix entry
#' @description Set a single source -> destination spillover coefficient in a
#'   spillover matrix by channel name, validating both channels against a panel
#'   table's `fcs_colname` column. Used to hand-tune a computed spillover matrix
#'   (e.g. CATALYST `computeSpillmat()` output) before applying compensation.
#' @param spill_matrix Numeric spillover matrix with channel names on both
#'   dimensions (rows = source channel, cols = destination channel).
#' @param panel_table Data frame with an `fcs_colname` column listing valid
#'   channel names; used only to validate `source_channel` / `dest_channel`.
#' @param source_channel Character, the emitting channel (row), e.g. "Gd157Di".
#' @param dest_channel Character, the receiving channel (column), e.g. "Gd156Di".
#' @param value Numeric spillover coefficient to set at
#'   `[source_channel, dest_channel]`.
#' @return The spillover matrix with the single entry updated.
#' @export
update_spill_matrix2 <- function(spill_matrix, panel_table, source_channel, dest_channel, value) {
  if (!source_channel %in% panel_table$fcs_colname | !dest_channel %in% panel_table$fcs_colname) {
    stop(paste0("Could not find channels ", source_channel, " and/or ", dest_channel, " in panel_table."))
  }
  spill_matrix[source_channel, dest_channel] <- value
  return(spill_matrix)
}
