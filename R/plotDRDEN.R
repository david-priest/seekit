# plotDRDEN.R — CATALYST-free rewrite of the CustomFunctionsSept25 function
# (my_functions.R). A faceted reduced-dimension plot (like plotDR2) with a 2D
# kernel-density (stat_density_2d) contour overlay per facet.
# Body verbatim from the legacy copy; only the CATALYST namespace shims
# (CATALYST:::.* internals, bare accessors) were rewritten to the package's
# .wl_* internals (R/wl_internals.R), ggrastr -> .wl_geom_point_rast, and the
# plotDR1_colorpal <<- global side-effect dropped (matches plotDR2).
# 2026-06 seekit migration of the CMV CyTOF pipeline.
plotDRDEN <- function (x, dr = NULL, textsize = 18, legendpointsize = 7, color_by = "condition", facet_by = NULL, hide_axis = F, alpha = 0.8, pointsize = 0.4,
                     ncol = NULL, assay = "exprs", scale = TRUE, random_order = FALSE, q = 0.01, dims = c(1, 2), k_pal = .wl_cluster_cols, 
                     a_pal = hcl.colors(10, "Viridis"), rast = FALSE, panel_spacing = 1, parse_color_by = FALSE, legend_title = NULL, 
                     bins = 10, density_alpha = 0.5, x_lim = NULL, y_lim = NULL) {
  stopifnot(is(x, "SingleCellExperiment"), .wl_check_assay(x, assay), length(reducedDims(x)) != 0, 
            is.logical(scale), length(scale) == 1, is.numeric(q), length(q) == 1, q >= 0, q < 0.5)
  
  .wl_check_pal(a_pal)
  .wl_check_cd_factor(x, facet_by)
  
  if (!is.null(ncol)) stopifnot(is.numeric(ncol), length(ncol) == 1, ncol%%1 == 0)
  
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
      nk <- nlevels(kids)
      if (length(k_pal) < nk) k_pal <- colorRampPalette(k_pal)(nk)
      
      # Save the color palette to a global variable
    } else {
      # Clustering is in colData
      kids <- factor(colData(x)[[color_by]])
    }
  }
  
  xy <- reducedDim(x, dr)[, dims]
  colnames(xy) <- c("x", "y")
  df <- data.frame(colData(x), xy, check.names = FALSE)
  
  if (!is.null(facet_by)) {
    stopifnot(facet_by %in% colnames(df))
    facet_levels <- unique(df[[facet_by]])
    
    
    # Create a list to store individual plots
    plot_list <- lapply(facet_levels, function(facet_level) {
      df_subset <- df[df[[facet_by]] == facet_level, ]
      
      # Shuffle rows of df_subset if random_order is TRUE
      if (random_order) {
        df_subset <- df_subset[sample(nrow(df_subset)), ]
      }
      
      # Apply the correct color palette
      if (!is.null(kids)) {
        df_subset[[color_by]] <- factor(df_subset[[color_by]], levels = levels(kids))
        color_scale <- scale_color_manual(values = k_pal)
      } else {
        color_scale <- scale_colour_gradientn(colors = a_pal)
      }
      
      # Use rasterized points if rast = TRUE
      if (rast) {
        point_layer <- .wl_geom_point_rast(size = pointsize, alpha = alpha, shape = 16, raster.dpi = 600)
      } else {
        point_layer <- geom_point(size = pointsize, alpha = alpha, shape = 16)
      }
      
      p <- ggplot(df_subset, aes(x = x, y = y, col = !!sym(color_by))) +
        point_layer +
        stat_density_2d(
          geom = "polygon",
          contour = TRUE,
          aes(fill = after_stat(level)),
          colour = "black",
          bins = bins,
          alpha = density_alpha
        ) +
        scale_fill_distiller(palette = "Blues", direction = 1) +
        color_scale +
        labs(
          title = facet_level  # Removed x and y labels
        ) +
        theme_minimal() +
        theme(
          text = element_text(size = textsize),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_blank(),  # Removed axis text
          axis.ticks = element_blank(),  # Removed axis ticks
          axis.title = element_blank(),  # Removed axis titles
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          panel.spacing = unit(panel_spacing, "lines"),
          aspect.ratio = 1
        )
      
      # Apply x and y axis limits if specified
      if (!is.null(x_lim)) {
        p <- p + xlim(x_lim)
      }
      if (!is.null(y_lim)) {
        p <- p + ylim(y_lim)
      }
      
      return(p)
    })
    
    # Combine plots using patchwork
    combined_plot <- patchwork::wrap_plots(plot_list, ncol = ncol) +
      patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom")
    
    return(combined_plot)
  } else {
    stop("facet_by must be specified to use this function.")
  }
}
