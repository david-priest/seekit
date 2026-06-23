# seekit

Shared R helpers for the Wing Lab — plotting, save helpers, SCE handling.

This package consolidates the helper `.R` files that used to be vendored as
per-project `R/<name>.R` copies across `~/My Drive/Wing Lab/`. One canonical
source, version-locked per project via `renv` for reproducibility.

## Install

In a project's R session:

```r
# one-time, per machine:
install.packages("remotes")

# pin to a specific commit/tag (reproducible):
remotes::install_github("david-priest/seekit@v0.0.0.9000")

# or floating to main (will pick up updates on next install):
remotes::install_github("david-priest/seekit")
```

For reproducibility, use `renv` in each analysis project:

```r
# in the project root:
renv::init()
renv::install("david-priest/seekit")
renv::snapshot()       # records the commit hash in renv.lock — commit it
```

When a helper gets fixed and you want this project to pick up the fix:

```r
renv::update("seekit")
renv::snapshot()       # commits the updated lock — old projects untouched
```

## Use

```r
library(seekit)
# every helper is exported under its function name:
f2(...)
plotDR1(...)
saveFig(...)
```

## Development workflow

1. Edit a helper here in `R/`.
2. Test locally: `devtools::load_all()` to source everything; `f2(...)` etc.
3. Commit + push. GitHub Actions runs `R CMD check` on every push and once a
   day to catch upstream CRAN breakage.
4. To cut a release: bump `Version:` in `DESCRIPTION`, tag the commit
   (`git tag v0.0.1; git push --tags`). Projects can pin to the tag.

## Importing helpers from the Wing Lab tree

The initial population of `R/` was done by
`tools/import_from_my_drive.R` — a re-runnable script that scans
`~/My Drive/Wing Lab/` for `**/R/*.R`, picks the newest-mtime copy of each
helper, and copies it in. Re-run to bring in newer in-the-wild versions
during the transition period (until every project has switched to importing
the package and stopped editing vendored copies):

```bash
Rscript tools/import_from_my_drive.R
```

The script lists the helpers it includes vs. skips (one-off project-specific
files like `gsea_*.R` are skipped — they don't belong in a shared package).

## Status

**Pre-roxygen.** Functions are exported via a catch-all `exportPattern` in
`NAMESPACE`. Adding `@param` / `@export` / `@return` roxygen comments per
function is a follow-up — once done, `devtools::document()` generates a clean
`NAMESPACE` and `man/` pages.

`R CMD check` will throw "no visible binding for global variable" warnings on
unqualified `dplyr::`/`ggplot2::` calls inside helpers. These are cosmetic
during the transition; cleanup is gradual.

## Related

- [`wing-lab-R-helpers`](../wing-lab-R-helpers/) — the older `wlrh` CLI that
  syncs vendored copies across projects. Use it for migrating helpers OUT of
  projects (replace `R/<name>.R` with `library(seekit)`) but
  `wlrh push` is now disabled by default to prevent surprise overwrites.
