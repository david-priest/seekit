# plot_spill
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 1562-1651). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

plot_spill <- function(spill_matrix, cytof_data, panel_table, n_cells = 1000, my_seed = 1234, ncol = 20, exclude_channels = NULL) {
  
  # Sub-sample cytof_data
  set.seed(my_seed)
  idx <- split(seq(ncol(cytof_data)), cytof_data$sample_id)
  idx <- lapply(idx, function(.) sample(., min(n_cells, length(.))))
  cytof_data <- cytof_data[, unlist(idx)]
  
  # Identify non-zero entries
  nonzero_entries <- which(spill_matrix != 0, arr.ind = TRUE)
  
  nonzout <<- nonzero_entries
  
  # Check if there are any non-zero entries
  if (nrow(nonzero_entries) == 0) {
    warning("No non-zero entries in the spill matrix. No plots to generate.")
    return(NULL)
  }
  
  # Extract numeric part from source channel names for ordering
  source_channels <- rownames(spill_matrix)[nonzero_entries[, "row"]]
  source_numbers <- as.numeric(gsub("[^0-9]", "", source_channels))
  
  # Order the entries by the numeric part of the source channel
  nonzero_entries <- nonzero_entries[order(source_numbers),]
  
  # Initialize an empty list to store the plots
  plot_list <- list()
  
  # Loop over non-zero entries and generate scatter plots
  for (i in 1:nrow(nonzero_entries)) {
    # Get the names of the source and destination channels
    source <- rownames(spill_matrix)[nonzero_entries[i, "row"]]
    dest <- colnames(spill_matrix)[nonzero_entries[i, "col"]]
    
    # If source and destination are same, skip to the next iteration
    if(source == dest) {
      next
    }
    
    # Exclude specified channels
    if (!is.null(exclude_channels) && (source %in% exclude_channels || dest %in% exclude_channels)) {
      next
    }
    
    # Get the corresponding antigens
    source_antigen <- panel_table$antigen[panel_table$fcs_colname == source]
    dest_antigen <- panel_table$antigen[panel_table$fcs_colname == dest]
    
    # Continue to the next iteration if any of the channels could not be mapped to antigens
    if (length(source_antigen) < 1 | length(dest_antigen) < 1) {
      message(paste0("Could not map channels ", source, " and/or ", dest, " to antigens. Skipping these channels."))
      next
    }
    
    # Check if antigens are valid row names in the cytof_data
    if (!(source_antigen %in% rownames(cytof_data)) | !(dest_antigen %in% rownames(cytof_data))) {
      message(paste0("Antigens ", source_antigen, " and/or ", dest_antigen, " are not valid in the CyTOF data. Skipping these antigens."))
      next
    }
    
    # Get spill value and format to 4 decimal places
    spill_value_raw <- spill_matrix[nonzero_entries[i, "row"], nonzero_entries[i, "col"]]
    spill_value <- sprintf("%.4f", spill_value_raw)
    
    # Create a title for the plot
    title <- paste0(source, " (", source_antigen, ") to ","\n", dest, " (", dest_antigen, ")",
                    "\nComp Value: ", spill_value)
    
    # Generate scatter plot and add it to the list
    p <- plotScatter(cytof_data, c(source_antigen, dest_antigen), zeros = T) + ggtitle(title)
    
    # Only add ggplot objects to the list
    if (inherits(p, "ggplot")) {
      plot_list[[i]] <- p
    }
  }
  
  # Remove NULL plots
  plot_list <- Filter(is.ggplot, plot_list)

  # If plot_list is not empty, create patchwork plot
  if(length(plot_list) > 0) {
    patchwork_plot <- wrap_plots(plot_list, ncol = ncol)
    return(patchwork_plot)
  } else {
    warning("No plots were generated.")
    return(NULL)
  }
}
