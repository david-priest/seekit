# wl_source_file.R
# -----------------------------------------------------------------------------
# Resolve the path of the source .Rmd/.qmd that called a save helper (f/f2/...),
# WITHOUT a hard dependency on rmdhelp. rmdhelp is GitHub-only (not on CRAN) and
# is not declared in DESCRIPTION, so a plain `library(rmdhelp)` blew up for users
# who installed seekit from the public repo. This helper degrades gracefully:
#   rmdhelp (if installed) -> knitr current input (during render) ->
#   rstudioapi active document (interactive) -> NA.
# Returns a single character path, or NA_character_ if it can't be determined.
# Callers must tolerate NA (fall back to a default folder name; skip the
# source-file copy).
.wl_source_file <- function() {
  tryCatch({
    if (requireNamespace("rmdhelp", quietly = TRUE)) {
      p <- tryCatch(rmdhelp::get_this_rmd_file(), error = function(e) NULL)
      if (length(p) == 1L && !is.na(p) && nzchar(p)) return(p)
    }
    if (requireNamespace("knitr", quietly = TRUE)) {
      ki <- tryCatch(knitr::current_input(dir = TRUE), error = function(e) NULL)
      if (length(ki) == 1L && !is.null(ki) && !is.na(ki) && nzchar(ki)) return(ki)
    }
    if (requireNamespace("rstudioapi", quietly = TRUE) &&
        isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
      rp <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                     error = function(e) NULL)
      if (length(rp) == 1L && !is.null(rp) && nzchar(rp)) return(rp)
    }
    NA_character_
  }, error = function(e) NA_character_)
}
