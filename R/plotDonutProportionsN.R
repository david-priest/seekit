# plotDonutProportionsN
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 5837-5930). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

plotDonutProportionsN <- function(sce, master_cluster_column, daughter_cluster_column, meta, sample_id_column, cols =  c(brewer.pal(9, "Paired"), "#9b59b6"), daughter_order = NULL, row_order = NULL, threshold = 10) {
  # Convert to data frame
  df <- as.data.frame(colData(sce))

  df[[master_cluster_column]] <- factor(df[[master_cluster_column]], levels = unique(df[[master_cluster_column]]))
  df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = unique(df[[daughter_cluster_column]])) # Ensure consistent factor levels
  df[[sample_id_column]] <- factor(df[[sample_id_column]], levels = unique(df[[sample_id_column]]))

  # Convert the daughter columns to factors with optional custom ordering
  if (!is.null(daughter_order)) {
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = daughter_order)
  } else {
    df[[daughter_cluster_column]] <- factor(df[[daughter_cluster_column]], levels = unique(df[[daughter_cluster_column]]))
  }

  # Calculate counts using dplyr to ensure only existing combinations are included
  count_df <- df %>%
    group_by_at(vars(sample_id_column, master_cluster_column, daughter_cluster_column)) %>%
    summarise(count = n(), .groups = "drop")

  # Ensure all combinations are present with zero counts where necessary
  complete_df <- count_df %>%
    complete(nesting(!!sym(sample_id_column), !!sym(master_cluster_column)), !!sym(daughter_cluster_column), fill = list(count = 0))

  # Calculate total counts for each master cluster within each sample
  master_cluster_counts <- complete_df %>%
    group_by_at(vars(sample_id_column, master_cluster_column)) %>%
    summarise(total_count = sum(count), .groups = "drop")

  # Label clusters that meet the threshold
  labeled_master_clusters <- master_cluster_counts %>%
    mutate(include_in_avg = total_count > threshold)

  # Merge the labels back to the complete_df
  labeled_df <- complete_df %>%
    left_join(labeled_master_clusters, by = c(sample_id_column, master_cluster_column))

  # Calculate proportions within each sample_id and master_cluster_column
  proportion_df <- labeled_df %>%
    group_by_at(vars(sample_id_column, master_cluster_column)) %>%
    mutate(proportion = count / sum(count) * 100) %>%
    ungroup()

  # Merge the cell counts with the proportion_df
  summary_df <- proportion_df %>%
    rename(cell_count = count)

  # Calculate the average proportions only for clusters that meet the threshold
  averaged_df <- summary_df %>%
    filter(include_in_avg) %>%
    group_by_at(vars(master_cluster_column, daughter_cluster_column)) %>%
    summarise(mean_proportion = mean(proportion), .groups = "drop") %>%
    ungroup()

  # Normalize the averaged proportions to ensure they sum to 100 within each master cluster
  normalized_df <- averaged_df %>%
    group_by_at(vars(master_cluster_column)) %>%
    mutate(normalized_proportion = mean_proportion / sum(mean_proportion) * 100) %>%
    ungroup()

  # Update the factor levels of master_cluster_column based on row_order
  if (!is.null(row_order)) {
    normalized_df[[master_cluster_column]] <- factor(normalized_df[[master_cluster_column]], levels = row_order)
  }

  # Create the plot
  p <- ggplot(normalized_df, aes(x = 2, y = normalized_proportion, fill = get(daughter_cluster_column))) +
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start = 0) +
    xlim(1, 2.5) +
    scale_fill_manual(values = cols, name = "Isotype") +
    theme_void() +
    facet_grid(rows = vars(get(master_cluster_column))) +
    theme(
      strip.text = element_text(size = 20), # Reduce strip text size
      legend.title = element_text(size = 20), # Reduce legend title size
      legend.text = element_text(size = 20), # Reduce legend text size
      panel.spacing = unit(0.1, "lines"), # Reduce space between panels
      plot.margin = margin(10, 10, 10, 10) # Reduce plot margins
    )

  # Return the final plot
  return(p)
}
