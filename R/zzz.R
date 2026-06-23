# Package load hook -----------------------------------------------------------
#
# A few helpers are drop-in clones of functions from heavier upstream packages
# (Seurat for DimPlot2/FeaturePlot2/VlnPlot2; shiny for cleanWorkspace) and reach
# the upstream functions by having their environment rebound to the upstream
# namespace -- this lets the function bodies call the upstream symbols unqualified
# without attaching the package or using `pkg:::internal`.
#
# That rebind MUST happen at load time, not at build/install time. Doing it as a
# top-level `environment(fn) <- asNamespace("Seurat")` evaluates the rebind while
# the lazy-load database is being prepared during `R CMD INSTALL`, which hard-
# fails on any machine (e.g. CI) that does not have the package installed:
#   Error in loadNamespace(name) : there is no package called 'Seurat'
#
# Performing the rebind in `.onLoad`, guarded by `requireNamespace`, means the
# package installs and `R CMD check`s cleanly without these (optional, Suggests)
# packages present, while still wiring up the upstream internals at runtime on a
# machine that does have them. The functions are only usable where their upstream
# is installed anyway, so the guard is behaviour-preserving.

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)

  rebind_to <- function(fns, pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) return(invisible(NULL))
    upstream <- asNamespace(pkg)
    for (fn in fns) {
      if (!exists(fn, envir = ns, inherits = FALSE)) next
      f <- get(fn, envir = ns)
      environment(f) <- upstream
      assign(fn, f, envir = ns)
    }
    invisible(NULL)
  }

  rebind_to(c("DimPlot2", "FeaturePlot2", "VlnPlot2", "ExIPlot2", "SingleExIPlot2"), "Seurat")
  rebind_to(c("cleanWorkspace"), "shiny")
}
