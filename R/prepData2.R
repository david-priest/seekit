#' @title Build a SingleCellExperiment from FCS data (CATALYST-free), with a Time option
#' @description Clean-room, CATALYST-free re-implementation of `CATALYST::prepData()`
#'   for the common mass-cytometry case (a `flowSet` or FCS directory + a `panel` and
#'   `md`). It builds the `counts` (and arcsinh `exprs`) assays, marker-class `rowData`,
#'   per-cell `colData` (sample_id + factors), and `experiment_info` metadata, moving
#'   non-mass channels (Time, Event_length, gaussian params, …) into `int_colData` and
#'   subsetting to mass channels — exactly like `prepData`.
#'
#'   The addition: **`keep_time`** surfaces the acquisition `Time` channel into the
#'   user-facing `colData` (column `time_col`), so you can colour/gate by it directly.
#'   Time is taken **raw** per event (CATALYST's cumulative cross-sample offset is
#'   intentionally omitted — it assumes sequential acquisitions and is wrong for
#'   debarcoded single-run data; add an offset yourself if you need monotonic time).
#'
#'   Simplifications vs CATALYST: `panel` and `md` are required (no `guessPanel`), and
#'   the FACS path / multi-FCS channel harmonisation are not reproduced.
#' @param x A `flowSet`, `flowFrame`, list of `flowFrame`s, or a character path to a
#'   directory of FCS files / vector of FCS paths.
#' @param panel,md Panel and metadata `data.frame`s (required).
#' @param features Channels to keep (default: the panel channels). `keep_time` also
#'   retains the Time channel regardless.
#' @param transform If `TRUE`, add an arcsinh `exprs` assay, Default: TRUE.
#' @param cofactor arcsinh cofactor, Default: 5.
#' @param panel_cols,md_cols Column-name maps, as in `CATALYST::prepData`.
#' @param by_time Order samples by `$BTIM` when `md` is `NULL`, Default: TRUE.
#' @param FACS If `TRUE`, keep all channels (no mass-channel subset), Default: FALSE.
#' @param keep_time Attach the raw acquisition Time to `colData`, Default: FALSE.
#' @param time_col Name of the new colData Time column, Default: 'Time'.
#' @param ... Passed to `read.FCS` when `x` is a path (e.g. `truncate_max_range`).
#' @return A `SingleCellExperiment`.
#' @export
#' @importFrom flowCore flowSet fsApply exprs identifier keyword read.flowSet isFCSfile
#' @importFrom SingleCellExperiment SingleCellExperiment int_colData int_colData<- int_metadata int_metadata<-
#' @importFrom SummarizedExperiment assay assay<- rowData colData
#' @importFrom S4Vectors DataFrame metadata
#' @export
prepData2 <- function(x, panel = NULL, md = NULL, features = NULL,
                      transform = TRUE, cofactor = 5,
                      panel_cols = list(channel = "fcs_colname", antigen = "antigen", class = "marker_class"),
                      md_cols = list(file = "file_name", id = "sample_id", factors = c("condition", "patient_id")),
                      by_time = TRUE, FACS = FALSE, keep_time = FALSE, time_col = "Time", ...) {

  if (is.null(panel) || is.null(md))
    stop("prepData2: `panel` and `md` are required (guessPanel is not reimplemented).", call. = FALSE)

  fs <- .wl_read_fs2(x, ...)

  panel <- data.frame(panel, check.names = FALSE, stringsAsFactors = FALSE)
  md    <- data.frame(md,    check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(c("channel", "antigen") %in% names(panel_cols),
            all(unlist(md_cols) %in% names(md)), c("file", "id", "factors") %in% names(md_cols))
  stopifnot(panel[[panel_cols$channel]] %in% flowCore::colnames(fs))

  # ---- channels to keep ----
  if (is.null(features)) features <- as.character(panel[[panel_cols$channel]])
  if (keep_time) {
    tch <- grep("time", flowCore::colnames(fs), ignore.case = TRUE, value = TRUE)
    if (length(tch) == 0)
      warning("prepData2: keep_time = TRUE but no Time channel in the flowSet.", call. = FALSE)
    features <- union(features, tch)
  }

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
  ant_panel[is.na(ant_panel)] <- as.character(panel[[panel_cols$channel]])[is.na(ant_panel)]

  # ---- subset to selected channels, then rename panel channels -> antigens ----
  # `hit[i]` = which panel row channel i is (0 for non-panel channels, e.g. Time).
  # Renaming per-channel via `hit` is order-robust (no reliance on fs/panel order).
  fs   <- fs[, features]
  chs0 <- flowCore::colnames(fs)                                  # original channel names (Cd106Di, …, Time)
  hit  <- match(chs0, panel[[panel_cols$channel]], nomatch = 0)
  ant  <- chs0
  ant[hit != 0] <- ant_panel[hit[hit != 0]]                       # panel channels -> antigen; others keep their name
  dup  <- table(ant)                                              # disambiguate duplicate antigens: a -> a.1, a.2
  for (a in names(dup)) if (dup[a] > 1) ant[ant == a] <- paste(a, seq_len(dup[a]), sep = ".")
  flowCore::colnames(fs) <- ant
  chs  <- flowCore::colnames(fs)

  # ---- counts matrix (channels x cells). NB: raw Time (no cumulative offset). ----
  es <- matrix(flowCore::fsApply(fs, flowCore::exprs), byrow = TRUE, nrow = length(chs), dimnames = list(chs, NULL))

  # ---- rowData: marker_class (non-panel channels e.g. Time -> "none") ----
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

  # ---- surface raw Time into user colData (before the mass-channel subset) ----
  if (keep_time) {
    ti <- grep("time", chs, ignore.case = TRUE)
    if (length(ti) == 1) sce[[time_col]] <- as.numeric(es[ti, ])
    else if (length(ti) > 1) warning("prepData2: multiple Time-like channels; not attaching.", call. = FALSE)
  }

  # ---- move non-mass channels to int_colData, subset to mass channels ----
  if (!FACS) {
    is_mass <- !is.na(.wl_ms_from_chs(chs0))
    icd <- S4Vectors::DataFrame(t(es[!is_mass, , drop = FALSE]), check.names = FALSE)  # cols already named by chs
    SingleCellExperiment::int_colData(sce) <- cbind(SingleCellExperiment::int_colData(sce), icd)
    sce <- sce[is_mass, ]
  }

  # ---- arcsinh transform ----
  if (transform) {
    SummarizedExperiment::assay(sce, "exprs", FALSE) <-
      asinh(sweep(SummarizedExperiment::assay(sce, "counts"), 1, cofactor, "/"))
    SingleCellExperiment::int_metadata(sce)$cofactor <- cofactor
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
