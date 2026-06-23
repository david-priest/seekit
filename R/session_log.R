#' @title Start a session log file
#' @description Tee R console output (stdout) and conditions
#'   (`message()`, `warning()`, `stop()`) into a log file while still
#'   displaying them in the IDE console. Designed for agent-assisted
#'   debugging: the agent reads the log file directly so the user
#'   doesn't have to copy-paste error messages.
#'
#' Implementation: `sink(con, type = "output", split = TRUE)` for
#' stdout (visible in both console and file). For conditions
#' (messages, warnings, errors), `globalCallingHandlers()` runs
#' lightweight handlers that copy `conditionMessage()` to the file
#' WITHOUT invoking muffle restarts, so the condition still
#' propagates normally to the IDE.
#'
#' @param path Log file path. Default: `"session.log"` in the current
#'   working directory.
#' @param append `TRUE` to append to an existing log (preserve prior
#'   session output); `FALSE` to overwrite. Default: `TRUE`.
#' @return Invisibly returns the absolute path to the log file.
#' @seealso `stop_session_log()`, `session_log_path()`.
#' @export
#' @examples
#' \dontrun{
#' # At the top of an .Rmd / .qmd:
#' seekit::start_session_log("session.log")
#'
#' # ...run chunks...
#'
#' # At end of session (optional -- closed automatically at R exit):
#' seekit::stop_session_log()
#' }
start_session_log <- function(path = "session.log", append = TRUE) {
  if (!is.null(.session_log_env$con)) {
    message(sprintf(
      "Session log already open at: %s. Call stop_session_log() first.",
      .session_log_env$path))
    return(invisible(.session_log_env$path))
  }
  abs_path <- normalizePath(path, mustWork = FALSE)
  dir.create(dirname(abs_path), recursive = TRUE, showWarnings = FALSE)

  # If `append = FALSE`, truncate the file first; subsequent writes
  # all happen in append mode (sink-via-connection + per-call
  # cat(file = path, append = TRUE) from condition handlers, both
  # targeting the same path).
  if (!append) file.create(abs_path, showWarnings = FALSE)

  # ONE connection for the stdout sink. Condition handlers don't
  # share this connection -- they open the same file fresh on each
  # call via cat(file = path, append = TRUE). Two long-lived
  # connections to the same file lose writes (the kernel/buffer
  # state on each handle can race); per-call append is reliable
  # and the volume of conditions is small.
  con_sink <- file(abs_path, open = "at")

  # Header (via con_sink -> appears in stdout too once the sink is
  # active; written before the sink call to land in the file with
  # the right ordering).
  cat(sprintf("\n\n# ── Session log %s ── #\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      file = con_sink)
  cat("# wd:    ", getwd(), "\n",
      "# R:     ", R.version.string, "\n",
      "# user:  ", Sys.info()[["user"]], "\n",
      "# host:  ", Sys.info()[["nodename"]], "\n",
      sep = "", file = con_sink)
  flush(con_sink)

  # stdout: file + console
  sink(con_sink, type = "output", split = TRUE)

  # Conditions: write to the path (open/close per call), no shared
  # connection. Console retains the normal condition propagation
  # because we do NOT invokeRestart() inside the handlers.
  globalCallingHandlers(
    message = .session_log_handler_msg,
    warning = .session_log_handler_warn,
    error   = .session_log_handler_err
  )

  .session_log_env$con  <- con_sink   # for closing
  .session_log_env$path <- abs_path   # handlers cat(file = path)

  message(sprintf("Session log started: %s", abs_path))
  invisible(abs_path)
}

#' @title Stop the active session log
#' @description Closes the file connection opened by
#'   `start_session_log()` and removes the condition handlers it
#'   installed. Safe to call when no log is active (no-op + message).
#' @return Invisibly returns the closed log path, or `NULL` if none.
#' @seealso `start_session_log()`.
#' @export
stop_session_log <- function() {
  if (is.null(.session_log_env$con)) {
    message("No active session log.")
    return(invisible(NULL))
  }

  # NOTE: we do NOT try to remove the globalCallingHandlers we
  # registered -- there is no documented way to unset them in R 4.x
  # (`globalCallingHandlers(message = NULL)` errors with "condition
  # handlers must be functions"). The handlers are no-ops when
  # `.session_log_env$con` is NULL, so leaving them installed is
  # harmless. They self-cleanup on R session exit.

  # Pop our stdout sink. Loop in case other sinks are stacked above.
  while (sink.number(type = "output") > 0) sink(NULL, type = "output")

  path <- .session_log_env$path
  cat(sprintf("\n# ── Session log closed %s ── #\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      file = .session_log_env$con)
  close(.session_log_env$con)
  .session_log_env$con  <- NULL
  .session_log_env$path <- NULL

  message(sprintf("Session log closed: %s", path))
  invisible(path)
}

#' @title Path of the active session log
#' @description Returns the absolute path of the currently active
#'   session log, or `NULL` if none.
#' @return Character path or `NULL`.
#' @seealso `start_session_log()`, `stop_session_log()`.
#' @export
session_log_path <- function() {
  if (is.null(.session_log_env$con)) NULL else .session_log_env$path
}

# ── Internal state + handlers ───────────────────────────────────────
# Hidden env so the log connection survives function scope and isn't
# stomped by interactive globalenv copies of these helpers.
.session_log_env <- new.env(parent = emptyenv())
.session_log_env$con  <- NULL   # sink target (split = TRUE)
.session_log_env$path <- NULL   # handlers cat(file = path, append = TRUE)

# Handlers open/close the file on each call. Two long-lived
# connections to one file lose writes due to per-handle buffer races;
# per-call append is reliable. Condition volume is low so the perf
# cost is negligible.
.session_log_handler_msg <- function(m) {
  if (!is.null(.session_log_env$path))
    cat(conditionMessage(m),
        file = .session_log_env$path, append = TRUE, sep = "")
}
.session_log_handler_warn <- function(w) {
  if (!is.null(.session_log_env$path))
    cat("\n[WARN] ", conditionMessage(w), "\n",
        file = .session_log_env$path, append = TRUE, sep = "")
}
.session_log_handler_err <- function(e) {
  if (!is.null(.session_log_env$path))
    cat("\n[ERROR] ", conditionMessage(e), "\n",
        file = .session_log_env$path, append = TRUE, sep = "")
}
