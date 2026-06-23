# plotDonutProportions
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 6211-6351). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

plotDonutProportions <- function(sce, master_cluster_column, daughter_cluster_column, meta, facet_column = NULL, sample_id_column, plot_order_df = NULL, daughter_order = NULL, facet_order_df = NULL, cell_threshold = 100) {
  # Convert to data frame
  df <- as.data.frame(colData(sce))

  # Convert the master, daughter, and sample_id columns to factors
  df[[master_cluster_column]] <- factor(df[[master_cluster_column]], levels = unique(df[[master_cluster_column]]))
  df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = unique(df[[daughter_cluster_column]])) # Ensure consistent factor levels
  df[[sample_id_column]] <- factor(df[[sample_id_column]], levels = unique(df[[sample_id_column]]))

  # Convert the facet column to a factor with optional custom ordering
  if (!is.null(facet_column)) {
    if (!is.null(facet_order_df)) {
      df[[facet_column]] <- factor(df[[facet_column]], levels = as.character(unlist(facet_order_df)))
    } else {
      df[[facet_column]] <- factor(df[[facet_column]], levels = unique(df[[facet_column]]))
    }
  }

  # Convert the daughter columns to factors with optional custom ordering
  if (!is.null(daughter_order)) {
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = daughter_order)
  } else {
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = unique(df[[daughter_cluster_column]]))
  }

  # Determine color palette based on unique levels of daughter clustering
  color_palette <- brewer.pal(min(length(levels(df[[daughter_cluster_column]])), 12), "Paired")

  # Calculate counts using dplyr to ensure only existing combinations are included
  if (!is.null(facet_column)) {
    count_df <- df %>%
      group_by_at(vars(sample_id_column, master_cluster_column, daughter_cluster_column, facet_column)) %>%
      summarise(count = n(), .groups = "drop")
  } else {
    count_df <- df %>%
      group_by_at(vars(sample_id_column, master_cluster_column, daughter_cluster_column)) %>%
      summarise(count = n(), .groups = "drop")
  }

  # Ensure all combinations are present with zero counts where necessary
  if (!is.null(facet_column)) {
    complete_df <- count_df %>%
      complete(nesting(!!sym(sample_id_column), !!sym(master_cluster_column), !!sym(facet_column)), !!sym(daughter_cluster_column), fill = list(count = 0))
  } else {
    complete_df <- count_df %>%
      complete(nesting(!!sym(sample_id_column), !!sym(master_cluster_column)), !!sym(daughter_cluster_column), fill = list(count = 0))
  }

  # Calculate proportions within each sample_id and master_cluster_column
  if (!is.null(facet_column)) {
    proportion_df <- complete_df %>%
      group_by_at(vars(sample_id_column, master_cluster_column, facet_column)) %>%
      mutate(proportion = count / sum(count) * 100) %>%
      ungroup()
  } else {
    proportion_df <- complete_df %>%
      group_by_at(vars(sample_id_column, master_cluster_column)) %>%
      mutate(proportion = count / sum(count) * 100) %>%
      ungroup()
  }

  # Merge the cell counts with the proportion_df
  summary_df <- proportion_df %>%
    rename(cell_count = count)

  # Output the dataframe to the global environment
  summary_out <<- summary_df
  
  # Calculate mean proportions and total cells
  if (!is.null(facet_column)) {
    averaged_df <- summary_df %>%
      group_by_at(vars(master_cluster_column, daughter_cluster_column, facet_column)) %>%
      summarise(mean_proportion = mean(proportion), total_cells = sum(cell_count), .groups = "drop") %>%
      ungroup()
  } else {
    averaged_df <- summary_df %>%
      group_by_at(vars(master_cluster_column, daughter_cluster_column)) %>%
      summarise(mean_proportion = mean(proportion), total_cells = sum(cell_count), .groups = "drop") %>%
      ungroup()
  }
  
  averaged_out <<- averaged_df
  
  # Filter out combinations of master_cluster_column and facet_column that have a sum of cells less than the threshold
  if (!is.null(facet_column)) {
    filtered_averaged_df <- averaged_df %>%
      group_by_at(vars(master_cluster_column, facet_column)) %>%
      filter(sum(total_cells) >= cell_threshold) %>%
      ungroup()
  } else {
    filtered_averaged_df <- averaged_df # Don't filter when not facetting
  }
  
  filtout <<- filtered_averaged_df
  
  # Normalize the averaged proportions to ensure they sum to 100 within each master cluster and facet
  if (!is.null(facet_column)) {
    normalized_df <- filtered_averaged_df %>%
      group_by_at(vars(master_cluster_column, facet_column)) %>%
      mutate(normalized_proportion = mean_proportion / sum(mean_proportion) * 100) %>%
      ungroup()
  } else {
    normalized_df <- filtered_averaged_df %>%
      group_by_at(vars(master_cluster_column)) %>%
      mutate(normalized_proportion = mean_proportion / sum(mean_proportion) * 100) %>%
      ungroup()
  }

  # Determine the order of master clusters for plotting
  plot_order <- if (!is.null(plot_order_df)) {
    as.character(unlist(plot_order_df))
  } else {
    levels(normalized_df[[master_cluster_column]])
  }

  # Create the plot
  p <- ggplot(normalized_df, aes(x = 2, y = normalized_proportion, fill = get(daughter_cluster_column))) +
    geom_bar(width = 1, stat = "identity",color = "black", size = 0.2) +
    coord_polar("y", start = 0) +
    xlim(1.5, 2.5) +
    scale_fill_manual(values = color_palette, name = "Isotype") +
    theme_void() +
    theme(
      strip.text.x = element_text(size = 20, angle = 90), # Rotate top strip text 90 degrees
      strip.text.y = element_text(size = 20), # Keep side strip text as is
      legend.title = element_text(size = 20), # Reduce legend title size
      legend.text = element_text(size = 20), # Reduce legend text size
      panel.spacing = unit(0.1, "lines"), # Reduce space between panels
      plot.margin = margin(10, 10, 10, 10) # Reduce plot margins
    )

  # Add faceting
  if (!is.null(facet_column)) {
    p <- p + facet_grid(rows = vars(get(master_cluster_column)), cols = vars(get(facet_column)))
  } else {
    p <- p + facet_wrap(vars(get(master_cluster_column)), ncol = 30) + theme(strip.text = element_text(size = 20))
  }

  # Return the final plot
  return(p)
}
