plotScatterGateGrid <- function (x, chs, color_by = NULL, facet_by = NULL, bins = 100,
                             assay = "exprs", label = c("target", "channel", "both"),
                             zeros = FALSE, gate_name = NULL, density_overlay = FALSE,
                             density_bins = 15, density_alpha = 0.3, use_points = FALSE,
                             point_pal = NULL, contour_pal = NULL, point_alpha = 0.2,
                             scales = "fixed", independent = "none",
                             facet_title = c("x", "y")) 
{
  if (!requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("Package \"ggnewscale\" needed for this function to work. Please install it.")
  }
  
  label <- match.arg(label)
  
  m <- rownames(x)
  c <- channels(x)
  i <- lapply(list(m, c), function(u) {
    i <- match(chs, u, nomatch = 0)
    if (all(i == 0)) NULL else i
  })
  i <- unlist(i)
  y <- x[i, , drop = FALSE]
  y <- assay(y, assay)
  nms <- switch(label, target = m, channel = c, both = ifelse(c == 
                                                                m, c, paste(c, m, sep = "-")))
  chs[i != 0] <- rownames(y) <- nms[i]
  if (isTRUE(color_by %in% names(.wl_cluster_codes(x)))) 
    x[[color_by]] <- .wl_cluster_ids(x, color_by)
  
  cd <- cbind(colData(x), int_colData(x))
  df <- data.frame(t(as.matrix(y)), cd, check.names = FALSE, 
                   stringsAsFactors = FALSE)
  cd_vars <- intersect(names(cd), names(df))
  
  if (length(chs) > 2) {
    if (!is.null(gate_name)) stop("Gating overlay is only supported for 2-channel scatter plots.")
    df <- reshape2::melt(df, id.vars = unique(c(chs[1], cd_vars)))
    facet <- "variable"
    ylab <- ylab(NULL)
    chs[2] <- "value"
  }
  else {
    facet <- NULL
    ylab <- NULL
  }
  
  xy <- sprintf("`%s`", chs)
  if (!zeros) 
    df <- df[rowSums(df[, chs[c(1, 2)]] == 0) == 0, ]
  
  facet_title <- match.arg(facet_title)
  facet_vars <- c(facet, facet_by)

  if (!is.null(facet_vars)) {
    wrap_formula <- as.formula(paste("~", paste(facet_vars, collapse = " + ")))

    # When facet_by has 2 variables, suppress the label for the unwanted direction
    if (!is.null(facet_by) && length(facet_by) == 2 && is.null(facet)) {
      keep_var   <- if (facet_title == "x") facet_by[1] else facet_by[2]
      drop_var   <- if (facet_title == "x") facet_by[2] else facet_by[1]
      lbl_list   <- setNames(list(label_value, function(x) rep("", length(x))),
                             c(keep_var, drop_var))
      lbl <- do.call(labeller, lbl_list)
    } else {
      lbl <- label_value
    }

    if (requireNamespace("ggh4x", quietly = TRUE) && independent != "none") {
      facet_layer <- ggh4x::facet_wrap2(wrap_formula, scales = scales, axes = "all", labeller = lbl)
    } else {
      facet_layer <- facet_wrap(wrap_formula, scales = scales, labeller = lbl)
    }
  } else {
    facet_layer <- NULL
  }
  
  if (is.null(color_by)) {
    col_var_points <- NULL
  } else {
    col_var_points <- sprintf("`%s`", color_by)
  }
  
  fill_var <- if (is.null(color_by)) "..ncount.." else NULL
  
  categorical_col_by <- !is.null(color_by) && !is.numeric(df[[color_by]])
  numeric_col_by <- !is.null(color_by) && is.numeric(df[[color_by]])
  
  if (categorical_col_by) {
    n_factors <- length(unique(df[[color_by]]))
    if (is.null(point_pal)) point_pal <- scales::hue_pal()(n_factors)
    if (is.null(contour_pal)) contour_pal <- scales::hue_pal(l = 40)(n_factors)
    
    if (is.null(names(point_pal))) names(point_pal) <- unique(df[[color_by]])
    if (is.null(names(contour_pal))) names(contour_pal) <- unique(df[[color_by]])
  } else if (numeric_col_by) {
    if (is.null(point_pal)) point_pal <- c("navy", rev(RColorBrewer::brewer.pal(11, "Spectral")))
    if (is.null(contour_pal)) contour_pal <- c("darkblue", rev(RColorBrewer::brewer.pal(11, "RdYlBu")))
  }
  
  p <- ggplot(df, aes_string(x = xy[1], y = xy[2], fill = fill_var)) +
    facet_layer + ylab + theme_bw() +
    theme(aspect.ratio = 1, panel.grid = element_blank(),
          axis.text = element_text(color = "black"), strip.background = element_rect(fill = "white"),
          legend.key.height = unit(0.8, "lines"))
  
  cat_guide_obj <- guide_legend(override.aes = list(alpha = 1, size = 3), order = 1)
  
  if (is.null(color_by)) {
    if (use_points) {
      p <- p + geom_point(alpha = point_alpha, size = 0.8, na.rm = TRUE, color = "black")
    } else {
      p <- p + geom_hex(bins = bins, na.rm = TRUE, show.legend = FALSE) +
        scale_fill_gradientn(trans = "sqrt", colors = c("navy", rev(RColorBrewer::brewer.pal(11, "Spectral"))))
    }
  } else {
    if (categorical_col_by) {
      p <- p + geom_point(aes_string(color = col_var_points), alpha = point_alpha, size = 0.8, na.rm = TRUE) +
        scale_color_manual(values = point_pal, guide = cat_guide_obj)
    } else if (numeric_col_by) {
      p <- p + geom_point(aes_string(color = col_var_points), alpha = point_alpha, size = 0.8, na.rm = TRUE) +
        scale_color_gradientn(colors = point_pal)
    }
  }
  
  if (density_overlay) {
    grp_cols_dens <- unique(c(facet_vars, color_by))
    if (length(grp_cols_dens) > 0) {
      df_dens <- df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols_dens))) |>
        dplyr::filter(dplyr::n() >= 5) |> 
        dplyr::filter(var(!!sym(chs[1]), na.rm = TRUE) > 0 & var(!!sym(chs[2]), na.rm = TRUE) > 0) |>
        dplyr::ungroup()
    } else {
      df_dens <- df
    }
    
    if (nrow(df_dens) > 0) {
      if (is.null(color_by)) {
        p <- p + stat_density_2d(data = df_dens, 
                                 inherit.aes = TRUE, 
                                 geom = "polygon",
                                 fill = "black", 
                                 alpha = density_alpha, 
                                 bins = density_bins, 
                                 na.rm = TRUE,
                                 contour_var = "ndensity")
      } else {
        p <- p + ggnewscale::new_scale_fill()
        
        if (categorical_col_by) {
          p <- p + stat_density_2d(data = df_dens,
                                   aes_string(fill = col_var_points), 
                                   inherit.aes = TRUE, 
                                   geom = "polygon",
                                   alpha = density_alpha, 
                                   bins = density_bins,
                                   na.rm = TRUE,
                                   contour_var = "ndensity") +
            scale_fill_manual(name = sprintf("%s (density)", color_by), values = contour_pal) 
          
        } else if (numeric_col_by) {
          p <- p + stat_density_2d(data = df_dens,
                                   aes_string(fill = col_var_points),
                                   inherit.aes = TRUE,
                                   geom = "polygon",
                                   alpha = density_alpha,
                                   bins = density_bins,
                                   na.rm = TRUE,
                                   contour_var = "ndensity") +
            scale_fill_gradientn(name = sprintf("%s (density)", color_by), colors = contour_pal)
        }
      }
    }
  }
  
  has_gate <- !is.null(gate_name) && (gate_name %in% colnames(df))
  if (has_gate) {
    poly_name <- paste0(gate_name, "_polygon")
    if (!poly_name %in% names(S4Vectors::metadata(x))) {
      stop(sprintf("Gating polygon coordinates not found in metadata under name: %s", poly_name))
    }
    
    poly_df <- S4Vectors::metadata(x)[[poly_name]]
    if (!identical(poly_df[1, ], poly_df[nrow(poly_df), ])) {
      poly_df <- rbind(poly_df, poly_df[1, ])
    }
    colnames(poly_df) <- chs[1:2]
    
    gate_vec <- df[[gate_name]]
    if (is.logical(gate_vec)) df$is_gated <- gate_vec else df$is_gated <- gate_vec == "Gated"
    
    grp_cols_pct <- unique(c(facet_vars, color_by))
    if (length(grp_cols_pct) > 0) {
      text_df <- df |>
        dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols_pct))) |>
        dplyr::summarise(pct = mean(is_gated, na.rm = TRUE) * 100, .groups = "drop")
    } else {
      text_df <- data.frame(pct = mean(df$is_gated, na.rm = TRUE) * 100)
    }
    
    text_df$label <- sprintf("%.1f%%", text_df$pct)
    text_df[[chs[1]]] <- -Inf
    text_df[[chs[2]]] <- Inf
    
    if (categorical_col_by) {
      if (!is.null(facet_vars)) {
        text_df <- text_df |>
          dplyr::group_by(dplyr::across(dplyr::all_of(facet_vars))) |>
          dplyr::mutate(vjust = dplyr::row_number() * 1.5) |>
          dplyr::ungroup()
      } else {
        text_df$vjust <- seq_len(nrow(text_df)) * 1.5
      }
    } else {
      text_df$vjust <- 1.5
    }
    
    p <- p + geom_polygon(data = poly_df, aes_string(x = xy[1], y = xy[2]), 
                          inherit.aes = FALSE, fill = NA, color = "black", size = 0.5)
    
    if (categorical_col_by) {
      p <- p + ggnewscale::new_scale_color()
      p <- p + geom_text(data = text_df, 
                         aes_string(x = xy[1], y = xy[2], label = "label", color = col_var_points, vjust = "vjust"), 
                         inherit.aes = FALSE, hjust = -0.1, fontface = "bold", show.legend = FALSE) +
        scale_color_manual(values = point_pal, guide = "none") 
    } else {
      p <- p + geom_text(data = text_df, 
                         aes_string(x = xy[1], y = xy[2], label = "label", vjust = "vjust"), 
                         inherit.aes = FALSE, color = "black", hjust = -0.1, fontface = "bold")
    }
  }
  
  return(p)
}

