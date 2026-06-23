# ei2.R — CATALYST-free rewrite of ei2 (per-sample count + metadata table).
# Body is David's own; only CATALYST's sample_ids() -> .wl_sample_ids().
# -----------------------------------------------------------------------------
ei2 <- function(sce, meta = c("patient_id")) {
  ns <- table(sample_id = .wl_sample_ids(sce))
  df <- as.data.frame(ns)
  m  <- match(df$sample_id, sce$sample_id)
  for (i in meta) df[[i]] <- sce[[i]][m]
  df
}
