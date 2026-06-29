# plotAbundanceStackedAIM.R — CATALYST-free rewrite promoted from dev/catalyst_quarantine/.
# Body verbatim from the project-vendored copy; only the CATALYST namespace
# shims (CATALYST:::.* internals, bare accessors, the asNamespace('CATALYST')
# hack) were rewritten to the package's .wl_* internals (R/wl_internals.R).
# 2026-06 seekit migration of the CMV CyTOF pipeline.
# plotAbundanceStackedAIM
# 
# Stacked-abundances plot across two mergings, one grouping variable (CI) (lines 1328-1463 of source .Rmd).
# Migrated from CMV CyTOF Figures David.Rmd as part of repository reorganisation
# (see CMV_paper_analysis.Rmd / CMV_extra_analyses.Rmd / CMV_code_quarantine.Rmd).


plotAbundanceStackedAIM <- function (x, 
                              k_abundances = "merging1",
                              k_split = "merging2",
                              sample_order = NULL, 
                              meta = c("sample_id", "patient_id", "condition"), 
                              group_by = "condition", 
                              shape_by = NULL, 
                              text_size = 16,
                              n_cols = 4, 
                              wrap_cols = 1,
                              rotang = 45, 
                              k_pal = .wl_cluster_cols, 
                              average_across_samples = TRUE,
                              plot_ci = FALSE) 
{
  # Extract cluster IDs from colData
  cluster_ids_abundances <- colData(x)[[k_abundances]]
  cluster_ids_split <- colData(x)[[k_split]]
  
  # Split data by k_split
  split_levels <- unique(cluster_ids_split)
  plots <- list()
  
  for (split_level in split_levels) {
    subset_indices <- which(cluster_ids_split == split_level)
    
    # Subset the SingleCellExperiment object
    subset_x <- x[, subset_indices]
    subset_cluster_ids <- cluster_ids_abundances[subset_indices]
    
    ns <- table(cluster_id = subset_cluster_ids, sample_id = colData(subset_x)$sample_id)
    fq <- prop.table(ns, 2) * 100
    df <- as.data.frame(fq)
    m <- match(df$sample_id, colData(subset_x)$sample_id)
    for (i in meta) df[[i]] <- colData(subset_x)[[i]][m]
    
    dfout <<- df
    
    if (average_across_samples) {
      # Calculate the average frequencies
      df_avg <- df %>%
        group_by(cluster_id, !!sym(group_by)) %>%
        summarise(Freq = mean(Freq)) %>%
        ungroup()
      
      # Calculate the 95% CI
      df_ci <- df %>%
        group_by(cluster_id, !!sym(group_by)) %>%
        summarise(CI_lower = mean(Freq) - 1.96 * (sd(Freq) / sqrt(n())),
                  CI_upper = mean(Freq) + 1.96 * (sd(Freq) / sqrt(n()))) %>%
        ungroup()
      
      # Join the CI back to the averages
      df <- left_join(df_avg, df_ci, by = c("cluster_id", group_by))
      
      # Calculate cumulative sum for error bar positioning
      df <- df %>%
        arrange(desc(cluster_id)) %>%
        group_by(!!sym(group_by)) %>%
        mutate(cumFreq = cumsum(Freq)) %>%
        ungroup()
    }
    
    dfoutavg <<- df
    
    if (!is.null(sample_order)) {
      df$sample_id <- factor(df$sample_id, ordered = TRUE, levels = sample_order)
    }
    
        # Save the data frame to the global environment before filtering
    assign(paste0("df_", split_level), df, envir = .GlobalEnv)
    
    # Filter out rows with NA in the group_by column
    df <- df %>%
      filter(!is.na(!!sym(group_by)))
    
    if (average_across_samples) {
      p <- ggplot(df, aes_string(x = group_by, y = "Freq", fill = "cluster_id")) + 
        labs(x = NULL, y = "Proportion [%]", title = paste(split_level)) + 
        theme_bw() + 
        theme(panel.grid = element_blank(), 
              strip.text = element_text(face = "bold"), 
              strip.background = element_rect(fill = NA, color = NA), 
              axis.text = element_text(color = "black"), 
              axis.ticks = element_line(color = "black"),
              text = element_text(size = text_size),
              axis.text.x = element_text(angle = rotang, hjust = 1, vjust = 1), 
              legend.key.height = unit(0.8, "lines")) +
        geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) + 
        scale_fill_manual("cluster_id", values = k_pal) + 
        scale_x_discrete(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(scale = 1)) + 
        theme(panel.border = element_blank(), panel.spacing.x = unit(1, "lines"))
      
      if (plot_ci) {
        p <- p + geom_errorbar(aes_string(ymin = "cumFreq - (Freq - CI_lower)", ymax = "cumFreq + (CI_upper - Freq)"), 
                               width = 0.5, position = position_dodge(width = 0.9))
      }
    } else {
      p <- ggplot(df, aes_string(x = "sample_id", y = "Freq", fill = "cluster_id")) + 
        labs(x = NULL, y = "Proportion [%]", title = paste(split_level)) + 
        theme_bw() + 
        theme(panel.grid = element_blank(), 
              strip.text = element_text(face = "bold"), 
              strip.background = element_rect(fill = NA, color = NA), 
              axis.text = element_text(color = "black"), 
              axis.ticks = element_line(color = "black"),
              text = element_text(size = text_size),
              axis.text.x = element_text(angle = rotang, hjust = 1, vjust = 1), 
              legend.key.height = unit(0.8, "lines")) +
        facet_wrap(group_by, scales = "free_x", ncol = n_cols) + 
        geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) + 
        scale_fill_manual("cluster_id", values = k_pal) + 
        scale_x_discrete(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(scale = 1)) + 
        theme(panel.border = element_blank(), panel.spacing.x = unit(1, "lines"))
    }
    
    p <- p + scale_fill_manual(values = k_pal)
    
    plots[[split_level]] <- p
  }
  
  combined_plot <- patchwork::wrap_plots(plots, ncol = wrap_cols)
  
  return(list(combined_plot = combined_plot, plots = plots))
}

