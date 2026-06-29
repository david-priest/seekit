# plotDRlite.R — CATALYST-free rewrite promoted from dev/catalyst_quarantine/.
# Body verbatim from the project-vendored copy; only the CATALYST namespace
# shims (CATALYST:::.* internals, bare accessors, the asNamespace('CATALYST')
# hack) were rewritten to the package's .wl_* internals (R/wl_internals.R).
# 2026-06 seekit migration of the CMV CyTOF pipeline.

plotDRlite <- function(
    x,
    dr                = NULL,
    textsize          = 18,
    legendpointsize   = 7,
    color_by          = "condition",
    facet_by          = NULL,
    hide_axis         = FALSE,
    border_width      = 1,
    ncol              = NULL,
    pointsize         = 0.4,
    assay             = "exprs",
    scale             = TRUE,
    random_order      = FALSE,
    q                 = 0.01,
    dims              = c(1, 2),
    alpha             = 0.8,
    k_pal             = .wl_cluster_cols,
    a_pal             = hcl.colors(10, "Viridis"),
    rast              = FALSE,
    raster.dpi        = 300,    # was hardcoded at 600; lower = faster render
    downsample        = NULL,   # integer: randomly subsample to this many cells after scaling
    panel_spacing     = 1,
    plot_order        = NULL,
    highlight_cluster = NULL
) {
    stopifnot(
        is(x, "SingleCellExperiment"),
        .wl_check_assay(x, assay),
        length(reducedDims(x)) != 0,
        is.logical(scale), length(scale) == 1,
        is.numeric(q), length(q) == 1, q >= 0, q < 0.5
    )

    .wl_check_pal(a_pal)
    .wl_check_cd_factor(x, facet_by)

    if (!is.null(ncol)) stopifnot(is.numeric(ncol), length(ncol) == 1, ncol %% 1 == 0)

    if (is.null(dr)) {
        dr <- reducedDimNames(x)[1]
    } else {
        stopifnot(is.character(dr), length(dr) == 1, dr %in% reducedDimNames(x))
    }

    stopifnot(is.numeric(dims), length(dims) == 2, dims %in% seq_len(ncol(reducedDim(x, dr))))

    if (!all(color_by %in% rownames(x))) {
        stopifnot(length(color_by) == 1)
        if (!color_by %in% names(colData(x))) {
            .wl_check_sce(x, TRUE)
            .wl_check_pal(k_pal)
            .wl_check_k(x, color_by)
            kids <- .wl_cluster_ids(x, color_by)
            nk   <- nlevels(kids)
            if (length(k_pal) < nk) k_pal <- colorRampPalette(k_pal)(nk)
            plotDR1_colorpal        <<- k_pal
            names(plotDR1_colorpal) <<- levels(kids)
        } else {
            kids <- NULL
        }
    }

    xy <- reducedDim(x, dr)[, dims]
    colnames(xy) <- c("x", "y")
    df <- data.frame(colData(x), xy, check.names = FALSE)

    if (all(color_by %in% rownames(x))) {
        color_by <- unique(color_by)

        # ---- Scale on ALL cells first (quantiles must reflect the full dataset) ----
        es <- as.matrix(assay(x, assay)[color_by, , drop = FALSE])
        if (scale) es <- .wl_scale_exprs(es, 1, q)

        # ---- Now filter NA coords and downsample --------------------------------
        valid <- !is.na(df$x) & !is.na(df$y)
        df    <- df[valid, ]
        es    <- es[, valid, drop = FALSE]

        if (!is.null(downsample) && nrow(df) > downsample) {
            set.seed(1)
            keep <- sort(sample(nrow(df), downsample))
            df   <- df[keep, ]
            es   <- es[, keep, drop = FALSE]
        }

        # ---- Build long format via row-indexing (no reshape2, no factor issues) --
        n_cells      <- nrow(df)
        df           <- df[rep(seq_len(n_cells), times = length(color_by)), , drop = FALSE]
        rownames(df) <- NULL
        df$value     <- c(t(es))   # t(es) = cells × markers; c() col-major = marker1 all cells, then marker2...
        df$variable  <- factor(rep(color_by, each = n_cells), levels = color_by)

        l        <- switch(assay, exprs = "expression", assay)
        l        <- paste0("scaled\n"[scale], l)
        scale    <- scale_colour_gradientn(l, colors = a_pal)
        thm      <- guide <- NULL
        color_by <- "value"
        facet    <- facet_wrap("variable", ncol = ncol)

    } else if (is.numeric(df[[color_by]])) {

        df <- df[!is.na(df$x) & !is.na(df$y), ]

        if (!is.null(downsample) && nrow(df) > downsample) {
            set.seed(1)
            df <- df[sort(sample(nrow(df), downsample)), ]
        }

        if (scale) {
            vs <- as.matrix(df[[color_by]])
            df[[color_by]] <- .wl_scale_exprs(vs, 2, q)
        }
        l        <- paste0("scaled\n"[scale], color_by)
        scale    <- scale_colour_gradientn(l, colors = a_pal)
        color_by <- sprintf("`%s`", color_by)
        facet    <- thm <- guide <- NULL

    } else {

        df <- df[!is.na(df$x) & !is.na(df$y), ]

        if (!is.null(kids)) kids <- kids[!is.na(xy[, 1]) & !is.na(xy[, 2])]

        if (!is.null(downsample) && nrow(df) > downsample) {
            set.seed(1)
            keep <- sort(sample(nrow(df), downsample))
            df   <- df[keep, ]
            if (!is.null(kids)) kids <- kids[keep]
        }

        facet <- NULL
        if (!is.null(kids)) {
            df[[color_by]] <- kids
            scale <- scale_color_manual(values = k_pal)
        } else {
            scale <- NULL
        }

        n     <- nlevels(droplevels(factor(df[[color_by]])))
        guide <- guides(col = guide_legend(
            ncol = ifelse(n > 12, 2, 1),
            override.aes = list(alpha = 1, size = legendpointsize)
        ))
        thm <- theme(legend.key.height = unit(0.8, "lines"), text = element_text(size = textsize))
    }

    if (dr %in% c("PCA", "MDS")) {
        asp <- coord_equal()
    } else {
        asp <- NULL
    }

    labs <- if (dr == "PCA") paste0("PC", dims) else paste(dr, "dim.", dims)

    # ---- point ordering (on already-filtered, already-downsampled df) -------
    if (!is.null(plot_order)) {
        df[[color_by]] <- factor(df[[color_by]], levels = plot_order)
        df <- df[order(df[[color_by]]), ]
    } else if (random_order) {
        set.seed(1)
        df <- df[sample(nrow(df)), ]
    }

    # ---- plot ---------------------------------------------------------------
    p <- ggplot(df, aes_string("x", "y", col = color_by))

    if (rast) {
        p <- p + geom_point_rast(size = pointsize, alpha = alpha, shape = 16, raster.dpi = raster.dpi)
    } else {
        p <- p + geom_point(size = pointsize, alpha = alpha, shape = 16)
    }

    if (!is.null(highlight_cluster)) {
        p <- p +
            geom_point(
                data = df[df[[color_by]] == highlight_cluster, ],
                aes(x = x, y = y),
                size = 2, alpha = 0.8, shape = 21, color = "black", fill = "yellow"
            ) +
            stat_density_2d(
                data = df[df[[color_by]] == highlight_cluster, ],
                aes(x = x, y = y, fill = after_stat(level)),
                geom = "polygon", color = "black", alpha = 0.3
            )
    }

    p <- p +
        labs(x = labs[1], y = labs[2]) +
        facet + scale + guide + asp +
        theme_minimal() +
        thm +
        theme(
            panel.grid.minor = element_blank(),
            strip.text       = element_text(face = "bold"),
            panel.grid.major = element_blank(),
            axis.text        = element_text(color = "black"),
            panel.spacing    = unit(panel_spacing, "lines"),
            aspect.ratio     = if (is.null(asp)) 1 else NULL
        ) +
        theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = border_width))

    if (hide_axis) {
        p <- p + theme(
            axis.text.x  = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank()
        )
    }

    if (is.null(facet_by)) return(p)

    if (is.null(facet)) {
        p + facet_wrap(facet_by, ncol = ncol)
    } else {
        if (nlevels(df$variable) == 1) {
            p + facet_wrap(facet_by, ncol = ncol) + ggtitle(levels(df$variable))
        } else {
            fs <- c("variable", facet_by)
            ns <- vapply(df[fs], nlevels, numeric(1))
            if (ns[2] > ns[1]) fs <- rev(fs)
            p + facet_grid(reformulate(fs[1], fs[2]))
        }
    }
}
