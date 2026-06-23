plotPbExprsCond <- function (x, k = "meta20", xaxis = "free_x", conditions = NULL, excluded_clusters = NULL,
                             features = "state", assay = "exprs", fun = c("median", "mean", "sum"), 
                             point_size = 1, squash = 0, clusters_order = NULL, textsize = 14, panel_spacing = 2, show_stats = TRUE, merging_col = F,
                             facet_by = c("antigen", "cluster_id"), color_by = "condition", mean_or_med = "median",
                             group_by = color_by, shape_by = NULL, size_by = FALSE, stat_size = 4,
                             geom = c("boxes","bar"), jitter = TRUE, ncol = NULL, group1keep = NULL, nudge = 0,
                             hide_x_labels = FALSE,
                             hide_ns = FALSE)   # TRUE hides non-significant brackets

{
  library(ggh4x) # Required for facet_wrap2
  library(rstatix)
  library(dplyr)
  library(ggplot2)
  
  fun <- match.arg(fun)
  geom <- match.arg(geom)
  facet_by <- match.arg(facet_by)
  stopifnot(is.logical(jitter), length(jitter) == 1)
  
  # Use the merging column from colData if provided
  if (!is.null(merging_col)) {
    cluster_ids <- factor(x[[k]])
  } else {
    # If no merging_col, proceed with normal clustering logic
    .wl_check_sce(x)
    k <- .wl_check_k(x, k)
    cluster_ids <- .wl_cluster_ids(x, k)
  }
  
  # Allow plotting from column data
  .wl_check_assay(x, assay)
  .wl_check_cd_factor(x, color_by)
  .wl_check_cd_factor(x, group_by)
  
  # Retrieve the features
  x <- x[.wl_get_features(x, features), ]
  x$cluster_id <- cluster_ids
  by <- c("cluster_id", "sample_id")
  ms <- .wl_agg(x, by, fun, assay)
  df <- melt(ms, varnames = c("antigen", by[length(by)]))
  
  if (length(by) == 2) 
    names(df)[ncol(df)] <- "cluster_id"
  x_var <- ifelse(facet_by == "antigen", group_by, "antigen")
  if (!is.null(df$cluster_id)) 
    df$cluster_id <- factor(df$cluster_id, levels(x$cluster_id))
  i <- match(df$sample_id, x$sample_id)
  j <- setdiff(names(colData(x)), c(names(df), "cluster_id"))
  df <- cbind(df, colData(x)[i, j, drop = FALSE])
  ncs <- table(as.list(colData(x)[by]))
  ncs <- rep(c(t(ncs)), each = nrow(x))
  if (size_by) {
    size_by <- "n_cells"
    df$n_cells <- ncs
  } else {
    size_by <- NULL
  }
  df <- df[ncs > 0, , drop = FALSE]
  
  if (!is.null(conditions)) {
    df <- df[df$condition %in% conditions, ]
    df$condition <- factor(df$condition, levels = conditions, ordered = TRUE)
  }
  
  # Add filtering for excluded_clusters (ok to exclude here because it's not proportions)
  if (!is.null(excluded_clusters)) {
    df <- df[!df$cluster_id %in% excluded_clusters, ]
  }
  
  if (!is.null(df$cluster_id) && !is.null(clusters_order)) {
    df$cluster_id <- factor(df$cluster_id, levels = clusters_order)
  }
  
  # Output the plotting data
  dfplotpbout <<- df
  
  # Calculate stats for all pairs using Wilcoxon with holm adjustment
  stat.test <- df %>%
    group_by(antigen, cluster_id) %>%
    wilcox_test(value ~ condition)  # For stats between conditions
  stat.test <- stat.test %>% add_y_position()
  
  # Keep only selected clusters in group1 if desired
  if (!is.null(group1keep)) {
    stat.test <- stat.test[stat.test$group1 %in% group1keep,]
    stat.test <- stat.test %>% add_y_position()
  }
  
  # Save the unfiltered stats to global
  pbCondStats <<- stat.test
  
  if (mean_or_med == "median") {
    fun = "median"
  } else {
    fun = "mean"
  }
  
  # Initialize ggplot
  p <- ggplot(df, aes_string(x = group_by, y =  "value")) + 
    facet_wrap2(facet_by, ncol = ncol, scales = xaxis, axes = "y") + 
    ylab(paste(features, "\n", fun, ifelse(assay == "exprs", "expression", assay))) +
    theme_bw() + 
    theme(
      panel.grid = element_blank(),
      text = element_text(size = textsize),
      strip.text = element_text(size = textsize),
      strip.background = element_rect(fill = NA, color = NA),
      strip.placement = "outside",  # Move facet labels outside
      legend.text = element_text(size = textsize),
      legend.title = element_text(size = textsize),
      panel.spacing = unit(panel_spacing, "lines"),
      axis.text = element_text(color = "black", size = textsize),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      #axis.text.x = element_blank(),   
      axis.title = element_text(size = textsize)
    )
  
  if (geom == "boxes") {
    p <- p + geom_boxplot(aes_string(fill = color_by), color = "black", width = 0.75, linewidth = 0.3, alpha = 0.9, outlier.color = NA, show.legend = T)
  } else if (geom == "bar") {
    p <- p + stat_summary(fun = fun, geom = "bar", aes_string(fill = color_by), color = "black", position = position_dodge(), alpha = 1)
  }
  
  # Define the position adjustment
  position_adjustment <- position_quasirandom(width = 0.2)
  
  # Add geom_quasirandom but no geom_segment in this function (one point per patient)
  p <- p + geom_quasirandom(width = 0.2, size = point_size, shape = 21, fill = "grey", color = "black", stroke = 0.5, show.legend = TRUE)
  
  if (hide_x_labels) {
    p <- p + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  }
  
  if (show_stats) {
    print("showing stats")
    p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.025, hide.ns = hide_ns, size = stat_size, bracket.nudge.y = nudge, step.increase = squash)
  }
  
  return(p)
}

