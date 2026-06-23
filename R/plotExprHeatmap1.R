# plotExprHeatmap1
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 282-437). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

plotExprHeatmap1 <- function (x, features = NULL, by = c("sample_id", "cluster_id", 
    "both"), k = "meta20", m = NULL, assay = "exprs", fun = c("median", 
    "mean", "sum"), scale = c("first", "last", "never"), q = 0.01, sub = FALSE,
    row_anno = TRUE, row_names = FALSE, col_anno = TRUE, row_clust = TRUE, col_clust = TRUE, 
    row_dend = TRUE, col_dend = TRUE, bars = FALSE, perc = FALSE, title = " ",
    bin_anno = FALSE, hm_pal = rev(brewer.pal(11, "RdYlBu")), 
    k_pal = .wl_cluster_cols, m_pal = k_pal, named_k_pal = NULL, distance = c("euclidean", 
        "maximum", "manhattan", "canberra", "binary", "minkowski"), 
    linkage = c("average", "ward.D", "single", "complete", "mcquitty", 
        "median", "centroid", "ward.D2"), plot_m_clusters = T) 
{
  
  # 
  if (sub == TRUE)
  {
    x <- subSCE(x, n_cells = 1000)
    title <- "WARNING: downsampled SCE used."
  }
  
    args <- as.list(environment())
    # .wl_check_args_plotExprHeatmap(args)
    distance <- match.arg(distance)
    linkage <- match.arg(linkage)
    scale <- match.arg(scale)
    fun <- match.arg(fun)
    by <- match.arg(by)
    x <- x[unique(.wl_get_features(x, features)), ]
    if (by != "sample_id") {
       .wl_check_k(x, k)
        x$cluster_id <- .wl_cluster_ids(x, k)
        #x$cluster_id <- x[[k]]
    }
    if (by == "both") 
        by <- c("cluster_id", "sample_id")
    .do_agg <- function() {
        z <- .wl_agg(x, by, fun, assay)
        if (length(by) == 1) 
            return(z)
        set_rownames(do.call("rbind", z), levels(x$cluster_id))
    }
    .do_scale <- function() {
        if (scale == "first") {
            z <- assay(x, assay)
            z <- .wl_scale_exprs(z, 1, q)
            assay(x, assay, FALSE) <- z
            return(x)
        }
        else .wl_scale_exprs(z, 1, q)
    }
    z <- switch(scale, first = {
        x <- .do_scale()
        .do_agg()
    }, last = {
        z <- .do_agg()
        .do_scale()
    }, never = {
        .do_agg()
    })
    if (length(by) == 1) 
        z <- t(z)
    if (scale != "never" && !(assay == "counts" && fun == "sum")) {
        qs <- round(quantile(z, c(0.01, 0.99)) * 5)/5
        lgd_aes <- list(at = seq(qs[1], qs[2], 0.2))
    }
    else lgd_aes <- list()
    lgd_aes$title_gp <- gpar(fontsize = 10, fontface = "bold", 
        lineheight = 0.8)
    if (!isFALSE(row_anno)) {
  if (plot_m_clusters) {
      left_anno <- switch(by[1], sample_id = .wl_anno_factors(x, levels(x$sample_id), row_anno, "row"), anno_clusters1(x, k, m, k_pal, m_pal, named_k_pal))
    } else {
      left_anno <- switch(by[1], sample_id = .wl_anno_factors(x, levels(x$sample_id), row_anno, "row"), anno_clusters1(x, k, NULL, k_pal, m_pal, named_k_pal))
    }
        catalyst_merge <<- merging_ids(x, k, m, k_pal, m_pal) # This will edit a global variable in the workspace with catalyst's choice for merging.
    }
    else left_anno <- NULL
    if (!isFALSE(col_anno) && length(by) == 2) {
        top_anno <- .wl_anno_factors(x, levels(x$sample_id), col_anno, 
            "colum")
    }
    else top_anno <- NULL
    if (bars) {
        right_anno <- anno_counts1(x[[by[1]]], perc)
    }
    else right_anno <- NULL
    if (bin_anno) {
        cell_fun <- function(j, i, x, y, ...) grid.text(gp = gpar(fontsize = 8), 
            sprintf("%.2f", z[i, j]), x, y)
    }
    else cell_fun <- NULL
    a <- ifelse(assay == "exprs", "expression", assay)
    f <- switch(fun, median = "med", fun)
    hm_title <- switch(scale, first = sprintf("%s %s\n%s", fun, 
        "scaled", a), last = sprintf("%s %s\n%s", "scaled", fun, 
        a), never = paste(fun, a, sep = "\n"))
    if (length(by) == 2) {
        col_title <- features
    }
    else if (length(features) == 1 && features %in% c("type", 
        "state")) {
        col_title <- paste0(features, "_markers")
    }
    else col_title <- ""
    
        if (!is.null(left_anno)) {
      heatmap_palette <<- left_anno@anno_list$cluster_id@color_mapping@colors
    }
  
    cn <- colnames(z) # Get the colunm names (names of each channel)
    
    p <- Heatmap(matrix = z, name = hm_title, 
        col = colorRamp2(seq(min(z), 
        max(z), l = n <- 100), 
        colorRampPalette(hm_pal)(n)), 
        column_title = col_title, column_title_side = ifelse(length(by) == 2, "top", "bottom"), 
        cell_fun = cell_fun, 
        cluster_rows = row_clust, 
        cluster_columns = col_clust, 
        show_row_dend = row_dend, 
        show_column_dend = col_dend, 
        clustering_distance_rows = distance, 
        clustering_method_rows = linkage, 
        clustering_distance_columns = distance, 
        clustering_method_columns = linkage, 
        show_row_names = row_names, # whether to show row names
        row_names_side = ifelse(by[1] == "cluster_id" || isFALSE(row_anno) && !row_dend || isFALSE(row_clust), "left", "right"), #which side for row names
        top_annotation = top_anno,  # Not sure what top annotation is
        left_annotation = left_anno, # Left annotation is the colours for clusters and merging (not row names)
        right_annotation = right_anno, # Right annotation is the bar plots and %
        rect_gp = gpar(col = "white"), # How to draw heatmap body colours
        heatmap_legend_param = lgd_aes,
        row_names_max_width = unit(13, "cm"),
        row_title = title,
        column_names_rot = 45,
        row_names_gp = gpar(fontsize = 12)) # Font size for row names (clusters)

  # Rotate the merging etc annotations too, but only if they exist
  if (!is.null(left_anno)) {
    p@left_annotation@anno_list[[1]]@name_param$rot = 45
    if (plot_m_clusters && length(p@left_annotation@anno_list) > 1) {
      p@left_annotation@anno_list[[2]]@name_param$rot = 45
    }
  }
    
    # Draw the heatmap
    ht <- draw(p)

    # Get the row labels in the order they appear on the clustered heatmap
    row_labels_after_clustering <- rownames(z)[row_order(ht)]
    row_labels_df <- data.frame(row_labels = row_labels_after_clustering)

    # Return both the heatmap object and the dataframe containing row labels
    list(heatmap = p, row_labels_df = row_labels_df)
}

