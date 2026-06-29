#' @title Build a SingleCellExperiment from FCS data (CATALYST-free), with Time/channel options
#' @description Clean-room, CATALYST-free re-implementation of `CATALYST::prepData()`
#'   for the common mass-cytometry case (a `flowSet` or FCS directory + a `panel` and
#'   `md`). It builds the `counts` (and `exprs`) assays, marker-class `rowData`,
#'   per-cell `colData` (sample_id + factors), and `experiment_info` metadata. By
#'   default non-mass channels (Time, Event_length, gaussian params, …) are moved into
#'   `int_colData` and the rows subset to mass channels — exactly like `prepData`.
#'
#'   Two additions over CATALYST:
#'   * **`keep_time`** surfaces the acquisition `Time` channel into the user-facing
#'     `colData` (column `time_col`), so you can colour/gate by it directly. Time is
#'     taken **raw** per event (CATALYST's cumulative cross-sample offset is omitted —
#'     it assumes sequential acquisitions and is wrong for debarcoded single-run data).
#'   * **`keep_channels`** retains chosen non-panel channels (Time, Event_length, the
#'     gaussian params, …) as **assay rows** (`marker_class = "none"`) instead of
#'     discarding them to `int_colData` — e.g. so a downstream gating app can gate on
#'     them. `keep_channels = "all"` keeps every non-mass channel. In `exprs`, kept
#'     gaussian params (Center/Width/Residual/Offset/Amplitude) are arcsinh-transformed
#'     like the metals (matching Cytobank / GateLabR gating space); only Time,
#'     Event_length, cell_length and file_number stay linear.
#'
#'   Simplifications vs CATALYST: `panel` and `md` are required (no `guessPanel`), and
#'   the multi-FCS channel harmonisation is not reproduced.
#' @param x A `flowSet`, `flowFrame`, list of `flowFrame`s, or a character path to a
#'   directory of FCS files / vector of FCS paths.
#' @param panel,md Panel and metadata `data.frame`s (required).
#' @param features Channels to keep (default: the panel channels). `keep_time` /
#'   `keep_channels` also retain their channels regardless.
#' @param transform If `TRUE`, add an `exprs` assay: `arcsinh(x/cofactor)` for all
#'   channels except the acquisition-level raw params (Time, Event_length, cell_length,
#'   file_number), which stay linear. Default: TRUE.
#' @param cofactor arcsinh cofactor (mass channels only), Default: 5.
#' @param panel_cols,md_cols Column-name maps, as in `CATALYST::prepData`.
#' @param by_time Order samples by `$BTIM` when `md` is `NULL`, Default: TRUE.
#' @param FACS If `TRUE`, keep all channels as rows (no mass-channel subset), Default: FALSE.
#' @param keep_time Attach the raw acquisition Time to `colData`, Default: FALSE.
#' @param time_col Name of the new colData Time column, Default: 'Time'.
#' @param keep_channels Non-panel channels to retain as assay rows (`marker_class
#'   "none"`, raw `exprs`): a character vector of channel names, or `"all"` for every
#'   non-mass channel. Default: NULL.
#' @param channel_names What to name the SCE rows (channels):
#'   * `"antigen"` (default) — rename panel channels to their antigen (`CD32`, …); the
#'     usual analysis form.
#'   * `"fcs"` — keep the raw FCS channel identifiers (`$PnN`: `Cd106Di`, `In115Di`, …).
#'     Use this for a pre-processing / debarcoding SCE (e.g. for GateLabR): the exported
#'     FCS keeps `$PnN`, so a later "real" import with a `$PnN`-keyed `panel` renames to
#'     antigens cleanly. No panel-based renaming is applied (marker_class still is).
#'   * `"desc"` — keep the raw FCS marker labels (`$PnS`: `106Cd_CD32`, `115In_CD45`, …;
#'     what GateLabR shows), falling back to `$PnN` where the desc is blank.
#' @param ... Passed to `read.FCS` when `x` is a path (e.g. `truncate_max_range`).
#' @return A `SingleCellExperiment`. `metadata()` also carries `pnn_to_channel` /
#'   `channel_to_pnn` (`$PnN` <-> marker/display-name maps) so gating apps (e.g.
#'   GateLabR) can resolve gating-ML metal-channel references and restore `$PnN` on
#'   FCS export.
#' @export
#' @importFrom flowCore flowSet fsApply exprs identifier keyword read.flowSet isFCSfile
#' @importFrom SingleCellExperiment SingleCellExperiment int_colData int_colData<- int_metadata int_metadata<-
#' @importFrom SummarizedExperiment assay assay<- rowData colData
#' @importFrom S4Vectors DataFrame metadata
prepData2 <- function(x, panel = NULL, md = NULL, features = NULL,
                      transform = TRUE, cofactor = 5,
                      panel_cols = list(channel = "fcs_colname", antigen = "antigen", class = "marker_class"),
                      md_cols = list(file = "file_name", id = "sample_id", factors = c("condition", "patient_id")),
                      by_time = TRUE, FACS = FALSE, keep_time = FALSE, time_col = "Time",
                      keep_channels = NULL, channel_names = c("antigen", "fcs", "desc"), ...) {

  channel_names <- match.arg(channel_names)

  if (is.null(panel) || is.null(md))
    stop("prepData2: `panel` and `md` are required (guessPanel is not reimplemented).", call. = FALSE)

  fs <- .wl_read_fs2(x, ...)

  panel <- data.frame(panel, check.names = FALSE, stringsAsFactors = FALSE)
  md    <- data.frame(md,    check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(c("channel", "antigen") %in% names(panel_cols),
            all(unlist(md_cols) %in% names(md)), c("file", "id", "factors") %in% names(md_cols))
  stopifnot(panel[[panel_cols$channel]] %in% flowCore::colnames(fs))

  all_chs   <- flowCore::colnames(fs)
  panel_chs <- as.character(panel[[panel_cols$channel]])

  # ---- resolve keep_channels -> kc (non-panel channels to retain as rows) ----
  if (FACS || identical(keep_channels, "all") || isTRUE(keep_channels)) {
    kc <- setdiff(all_chs, panel_chs)                       # every non-panel channel
  } else if (!is.null(keep_channels)) {
    kc   <- intersect(as.character(keep_channels), all_chs)
    miss <- setdiff(as.character(keep_channels), all_chs)
    if (length(miss))
      warning("prepData2: keep_channels not in flowSet: ", paste(miss, collapse = ", "), call. = FALSE)
  } else kc <- character(0)

  # ---- channels to read from fs ----
  if (is.null(features)) features <- panel_chs
  if (keep_time) {
    tch <- grep("time", all_chs, ignore.case = TRUE, value = TRUE)
    if (length(tch) == 0)
      warning("prepData2: keep_time = TRUE but no Time channel in the flowSet.", call. = FALSE)
    features <- union(features, tch)
  }
  features <- union(features, kc)

  # ---- match md files to flowSet, reorder fs ----
  ids0 <- md[[md_cols$file]]
  ids1 <- flowCore::fsApply(fs, flowCore::identifier)
  if (!all(ids1 %in% ids0))
    stop("prepData2: couldn't match flowSet/FCS filenames to md[[md_cols$file]].", call. = FALSE)
  fs <- fs[match(ids0, ids1)]

  # ---- md -> per-sample factor table, id renamed to sample_id ----
  k  <- c(md_cols$id, md_cols$factors)
  md <- md[, k, drop = FALSE]
  md[] <- lapply(md, factor)
  names(md)[1] <- "sample_id"

  # ---- antigen names (fall back to channel where the panel antigen is missing) ----
  ant_panel <- as.character(panel[[panel_cols$antigen]])
  ant_panel[is.na(ant_panel)] <- panel_chs[is.na(ant_panel)]

  # ---- subset to selected channels, then set row names per `channel_names` ----
  # `hit[i]` = which panel row channel i is (0 for non-panel channels, e.g. Time).
  # Renaming per-channel via `hit` is order-robust (no reliance on fs/panel order).
  fs   <- fs[, features]
  chs0 <- flowCore::colnames(fs)                                  # raw FCS channel names ($PnN: Cd106Di, …, Time)
  hit  <- match(chs0, panel_chs, nomatch = 0)
  ant  <- switch(channel_names,
    # panel antigens (CD32, …); non-panel channels keep their $PnN
    antigen = { a <- chs0; a[hit != 0] <- ant_panel[hit[hit != 0]]; a },
    # raw FCS channel identifiers ($PnN: Cd106Di, In115Di, …) — round-trips with a $PnN-keyed panel
    fcs     = chs0,
    # raw FCS marker labels ($PnS: 106Cd_CD32, 115In_CD45, …; GateLabR's display), $PnN where desc is blank
    desc    = {
      pp <- flowCore::pData(flowCore::parameters(fs[[1]]))
      dm <- stats::setNames(trimws(as.character(pp$desc)), as.character(pp$name))
      d  <- unname(dm[chs0]); na <- is.na(d) | d == ""; d[na] <- chs0[na]; d
    })
  dup  <- table(ant)                                              # disambiguate duplicate names: a -> a.1, a.2
  for (a in names(dup)) if (dup[a] > 1) ant[ant == a] <- paste(a, seq_len(dup[a]), sep = ".")
  flowCore::colnames(fs) <- ant
  chs  <- flowCore::colnames(fs)

  # ---- counts matrix (channels x cells). NB: raw Time (no cumulative offset). ----
  es <- matrix(flowCore::fsApply(fs, flowCore::exprs), byrow = TRUE, nrow = length(chs), dimnames = list(chs, NULL))

  # ---- rowData: marker_class (non-panel channels e.g. Time, gaussian -> "none") ----
  mc_lvls <- c("type", "state", "none")
  mcs_chr <- rep("none", length(chs))
  if (!is.null(panel_cols$class) && !is.null(panel[[panel_cols$class]])) {
    pc <- as.character(panel[[panel_cols$class]])
    if (any(!pc %in% mc_lvls))
      stop("prepData2: invalid marker classes; valid are 'type', 'state', 'none'.", call. = FALSE)
    mcs_chr[hit != 0] <- pc[hit[hit != 0]]
  }
  mcs <- factor(mcs_chr, levels = mc_lvls)
  rd <- S4Vectors::DataFrame(row.names = chs, channel_name = chs0, marker_name = chs, marker_class = mcs)

  # ---- colData: expand per-sample md to per-cell ----
  md$n_cells <- as.numeric(flowCore::fsApply(fs, nrow))
  kk <- setdiff(names(md), "n_cells")
  cd <- S4Vectors::DataFrame(lapply(md[kk], function(u)
    factor(as.character(rep(u, md$n_cells)), levels = levels(u))), row.names = NULL)

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = es), rowData = rd, colData = cd,
    metadata = list(experiment_info = .wl_build_ei(cd)))

  # ---- surface raw Time into user colData (before any row subset) ----
  if (keep_time) {
    ti <- grep("time", chs, ignore.case = TRUE)
    if (length(ti) == 1) sce[[time_col]] <- as.numeric(es[ti, ])
    else if (length(ti) > 1) warning("prepData2: multiple Time-like channels; not attaching.", call. = FALSE)
  }

  # ---- keep mass channels (+ any non-mass requested via keep_channels / FACS) as
  #      rows; move the remaining non-mass channels to int_colData ----
  is_mass  <- !is.na(.wl_ms_from_chs(chs0))
  keep_row <- is_mass | (chs0 %in% kc)
  move     <- !keep_row
  # Don't duplicate Time into int_colData when keep_time already surfaced it to colData —
  # else cbind(colData, int_colData) (e.g. inside CATALYST::plotScatter) has two `Time`
  # columns and newer ggplot2 errors ("data must be uniquely named ...").
  if (keep_time) move[grep("time", chs, ignore.case = TRUE)] <- FALSE
  if (any(move)) {
    icd <- S4Vectors::DataFrame(t(es[move, , drop = FALSE]), check.names = FALSE)  # cols named by chs
    SingleCellExperiment::int_colData(sce) <- cbind(SingleCellExperiment::int_colData(sce), icd)
  }
  sce <- sce[keep_row, ]

  # ---- transform: match Cytobank / GateLabR CyTOF display space — arcsinh(x/cofactor)
  #      for ALL channels (metal AND gaussian params: Center/Width/Residual/Offset/
  #      Amplitude) EXCEPT the acquisition-level raw params (Time, Event_length,
  #      cell_length, file_number), which stay linear. This keeps kept gaussian channels
  #      on the same scale gating-ML gates were drawn in. ----
  if (transform) {
    cnts <- SummarizedExperiment::assay(sce, "counts")
    ex   <- asinh(cnts / cofactor)
    raw  <- tolower(SummarizedExperiment::rowData(sce)$channel_name) %in%
            c("time", "event_length", "cell_length", "file_number")
    if (any(raw)) ex[raw, ] <- cnts[raw, ]
    SummarizedExperiment::assay(sce, "exprs", FALSE) <- ex
    SingleCellExperiment::int_metadata(sce)$cofactor <- cofactor
  }

  # ---- $PnN <-> display (rowname) maps, for gating apps (e.g. GateLabR): resolve
  #      gating-ML metal-channel refs (by $PnN) to the SCE rows, and restore $PnN on FCS
  #      export. Harmless in any channel_names mode. ----
  pnn  <- as.character(SummarizedExperiment::rowData(sce)$channel_name)   # $PnN
  disp <- rownames(sce)                                                   # marker/display name
  S4Vectors::metadata(sce)$pnn_to_channel <- as.list(stats::setNames(disp, pnn))  # $PnN -> display
  S4Vectors::metadata(sce)$channel_to_pnn <- as.list(stats::setNames(pnn, disp))  # display -> $PnN

  # ---- declare the transform so a gating app doesn't re-detect/re-scale (CyTOF arcsinh) ----
  if (transform) {
    S4Vectors::metadata(sce)$instrument_type <- "cytof"
    S4Vectors::metadata(sce)$transform_type  <- "arcsinh"
    S4Vectors::metadata(sce)$cofactor        <- cofactor
  }
  sce
}

# --- file-local helpers (CATALYST-free) -------------------------------------

# mass number from a channel name ("Cd106Di" -> 106; "Time" -> NA)
.wl_ms_from_chs <- function(chs) suppressWarnings(as.numeric(gsub("[[:punct:][:alpha:]]", "", chs)))

# read x into a flowSet (flowSet | flowFrame | list<flowFrame> | path(s))
.wl_read_fs2 <- function(x, ...) {
  if (methods::is(x, "flowSet")) return(x)
  if (methods::is(x, "flowFrame")) return(flowCore::flowSet(x))
  if (is.list(x) && all(vapply(x, methods::is, logical(1), "flowFrame")))
    return(flowCore::flowSet(x))
  if (is.character(x)) {
    fcs <- if (length(x) == 1 && dir.exists(x))
      list.files(x, pattern = "[.]fcs$", full.names = TRUE, ignore.case = TRUE) else x
    if (length(fcs) == 0) stop("prepData2: no FCS files found.", call. = FALSE)
    args <- list(...)
    for (. in c("transformation", "truncate_max_range"))
      if (is.null(args[[.]])) args[[.]] <- FALSE
    return(do.call(flowCore::read.flowSet, c(list(fcs), args)))
  }
  stop("prepData2: `x` must be a flowSet / flowFrame / list of flowFrames / FCS path(s).", call. = FALSE)
}

# experiment_info table (reimplements CATALYST:::.get_ei)
.wl_build_ei <- function(cd) {
  cd  <- as.data.frame(cd)
  ids <- levels(droplevels(factor(cd$sample_id)))
  j <- setdiff(names(cd), "cluster_id")
  j <- j[vapply(j, function(.) !is.numeric(cd[[.]]), logical(1))]
  j <- j[vapply(j, function(.) {
    ns <- table(cd$sample_id, cd[[.]]); all(rowSums(ns != 0) == 1)
  }, logical(1))]
  i   <- match(ids, cd$sample_id)
  ncs <- table(cd$sample_id)
  data.frame(cd[i, j, drop = FALSE], row.names = NULL, n_cells = as.integer(ncs))
}
