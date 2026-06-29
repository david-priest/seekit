# plotAbundanceStacked2.R — CATALYST-free rewrite promoted from dev/catalyst_quarantine/.
# Body verbatim from the project-vendored copy; only the CATALYST namespace
# shims (CATALYST:::.* internals, bare accessors, the asNamespace('CATALYST')
# hack) were rewritten to the package's .wl_* internals (R/wl_internals.R).
# 2026-06 seekit migration of the CMV CyTOF pipeline.
# plotAbundanceStacked2
# 
# Stacked-abundances plot across two mergings split by group_by and AIM_condition (lines 1470-1606 of source .Rmd).
# Migrated from CMV CyTOF Figures David.Rmd as part of repository reorganisation
# (see CMV_paper_analysis.Rmd / CMV_extra_analyses.Rmd / CMV_code_quarantine.Rmd).


plotAbundanceStacked2 <- function (x, 
                              k_abundances = "merging1",
                              k_split = "merging2",
                              sample_order = NULL, 
                              meta = c("sample_id", "patient_id", "condition", "AIM_cond"), 
                              group_by = "condition", 
                              shape_by = NULL, 
                              n_cols = 4, 
                              wrap_cols = 1,
                              rotang = 45, 
                              k_pal = .wl_cluster_cols, 
                              average_across_samples = TRUE,
                              text_size = 16,
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
      # Calculate the average frequencies within each AIM_cond
      df_avg <- df %>%
        group_by(cluster_id, !!sym(group_by), AIM_cond) %>%
        summarise(Freq = mean(Freq)) %>%
        ungroup()
      
      # Calculate the 95% CI within each AIM_cond
      df_ci <- df %>%
        group_by(cluster_id, !!sym(group_by), AIM_cond) %>%
        summarise(CI_lower = mean(Freq) - 1.96 * (sd(Freq) / sqrt(n())),
                  CI_upper = mean(Freq) + 1.96 * (sd(Freq) / sqrt(n()))) %>%
        ungroup()
      
      # Join the CI back to the averages
      df <- left_join(df_avg, df_ci, by = c("cluster_id", group_by, "AIM_cond"))
      
      # Calculate cumulative sum for error bar positioning
      df <- df %>%
        arrange(desc(cluster_id)) %>%
        group_by(!!sym(group_by), AIM_cond) %>%
        mutate(cumFreq = cumsum(Freq)) %>%
        ungroup()
    }
    
    dfoutavg <<- df
    
    if (!is.null(sample_order)) {
      df$sample_id <- factor(df$sample_id, ordered = TRUE, levels = sample_order)
    }
    
        # Filter out rows with NA in the group_by column
    df <- df %>%
      filter(!is.na(!!sym(group_by)))
    
    if (average_across_samples) {
      p <- ggplot(df, aes_string(x = group_by, y = "Freq", fill = "cluster_id")) + 
        facet_wrap(~AIM_cond, scales = "free_x", ncol = n_cols) +
        labs(x = NULL, y = "Proportion [%]", title = paste(split_level)) + 
        theme_bw() + 
        theme(panel.grid = element_blank(), 
              strip.text = element_text(face = "bold"), 
              strip.background = element_rect(fill = NA, color = NA), 
              axis.text = element_text(color = "black"), 
              axis.ticks = element_line(color = "black"),
              axis.text.x = element_text(angle = rotang, hjust = 1, vjust = 1), 
              text = element_text(size = text_size),
              legend.key.height = unit(0.8, "lines")) +
        geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) + 
        scale_fill_manual("cluster_id", values = k_pal) + 
        scale_x_discrete(expand = c(0, 0)) + 
        scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(scale = 1)) + 
        theme(panel.border = element_blank(), panel.spacing.x = unit(1, "lines"))
      
      if (plot_ci) {
        p <- p + geom_errorbar(aes_string(ymin = "cumFreq - (Freq - CI_lower)", ymax = "cumFreq + (CI_upper - Freq)"), 
                               width = 0.4, position = position_dodge(width = 0.9))
      }
    } else {
      df <- df %>%
        mutate(group_AIM = paste(!!sym(group_by), AIM_cond, sep = "_"))
      
      p <- ggplot(df, aes_string(x = "patient_id", y = "Freq", fill = "cluster_id")) + 
        facet_wrap(~group_AIM, scales = "free_x", ncol = n_cols) +
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
    }
    
    p <- p + scale_fill_manual(values = k_pal)
    
    plots[[split_level]] <- p
  }
  
  combined_plot <- patchwork::wrap_plots(plots, ncol = wrap_cols)
  
  return(list(combined_plot = combined_plot, plots = plots))
}

