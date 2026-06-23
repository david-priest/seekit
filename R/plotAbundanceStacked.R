# plotAbundanceStacked.R — CATALYST-free rewrite (.cluster_cols -> .wl_cluster_cols, namespace hack removed).
# plotAbundanceStacked
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 3031-3149). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

plotAbundanceStacked <- function (x, 
                                  k_abundances = "merging1",
                                  sample_order = NULL, 
                                  meta = c("sample_id", "patient_id", "condition"), 
                                  facet_by = NULL,
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
  
  ns <- table(cluster_id = cluster_ids_abundances, sample_id = colData(x)$sample_id)
  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)
  m <- match(df$sample_id, colData(x)$sample_id)
  for (i in meta) df[[i]] <- colData(x)[[i]][m]

  if (average_across_samples) {
    # Average (and 95% CI) grouped by the FULL facet unit: cluster_id x group_by x
    # facet_by. The earlier version grouped the mean by (cluster_id, group_by) only
    # and re-attached facet_by via an ei2() left_join, which DUPLICATES rows when
    # facet_by isn't nested in group_by. Grouping directly is correct and removes
    # the join. (Reconciled with the nBass_helpers local copy + kept the CI; the
    # two debug `dfout <<-`/`dfoutavg <<-` global assignments were removed.)
    grp <- c("cluster_id", group_by, if (!is.null(facet_by)) facet_by)
    df_avg <- df %>%
      group_by(across(all_of(grp))) %>%
      summarise(Freq = mean(Freq), .groups = "drop")
    df_ci <- df %>%
      group_by(across(all_of(grp))) %>%
      summarise(CI_lower = mean(Freq) - 1.96 * (sd(Freq) / sqrt(n())),
                CI_upper = mean(Freq) + 1.96 * (sd(Freq) / sqrt(n())),
                .groups = "drop")
    df2 <- left_join(df_avg, df_ci, by = grp)

    # Cumulative sum (for error-bar positioning) within each bar = group_by x facet_by.
    df2 <- df2 %>%
      arrange(desc(cluster_id)) %>%
      group_by(across(all_of(c(group_by, if (!is.null(facet_by)) facet_by)))) %>%
      mutate(cumFreq = cumsum(Freq)) %>%
      ungroup()
    df <- df2
  }

  if (!is.null(sample_order)) {
    df$sample_id <- factor(df$sample_id, ordered = TRUE, levels = sample_order)
  }
  
  if (average_across_samples) {
    p <- ggplot(df, aes_string(x = group_by, y = "Freq", fill = "cluster_id")) + 
      labs(x = NULL, y = "Proportion [%]") + 
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            strip.text = element_text(face = "bold"), 
            strip.background = element_rect(fill = "grey88", color = NA), 
            axis.text = element_text(color = "black"), 
            text = element_text(size = text_size),
            axis.text.x = element_text(angle = rotang, hjust = 1, vjust = 1), 
            legend.key.height = unit(0.8, "lines"),
            panel.spacing.x = unit(0.4, "lines")) +  # Adjust spacing between bars
      geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) + 
      scale_fill_manual("cluster_id", values = k_pal) + 
      scale_x_discrete(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(scale = 1)) + 
      theme(panel.border = element_blank(), panel.spacing.x = unit(1, "lines"))
    
    if (!is.null(facet_by)) {
      #p <- p + facet_wrap(facet_by, scales = "free_x", ncol = wrap_cols)
      #p <- p + facet_grid2(cols = vars(!!sym(facet_by)), scales = "free_x", space = "free")
      p <- p + facet_grid2(cols = vars(!!sym(facet_by)), scales = "free_x", space = "free") +
              theme(strip.placement = "outside")
    }
    
    if (plot_ci) {
      p <- p + geom_errorbar(aes_string(ymin = "cumFreq - (Freq - CI_lower)", ymax = "cumFreq + (CI_upper - Freq)"), 
                             width = 0.5, position = position_dodge(width = 0.9))
    }
  } else {
    p <- ggplot(df, aes_string(x = "sample_id", y = "Freq", fill = "cluster_id")) + 
      labs(x = NULL, y = "Proportion [%]") + 
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            strip.text = element_text(face = "bold"), 
            strip.background = element_rect(fill = "grey88", color = NA), 
            axis.text = element_text(color = "black"), 
            text = element_text(size = text_size),
            axis.text.x = element_text(angle = rotang, hjust = 1, vjust = 1), 
            legend.key.height = unit(0.8, "lines"),
            panel.spacing.x = unit(0.1, "lines")) +  # Adjust spacing between bars
      facet_wrap(facet_by, scales = "free_x", ncol = n_cols) + 
      geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) + 
      scale_fill_manual("cluster_id", values = k_pal) + 
      scale_x_discrete(expand = c(0, 0)) + 
      scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(scale = 1)) + 
      theme(panel.border = element_blank(), panel.spacing.x = unit(1, "lines"))
  }
  
  p <- p + scale_fill_manual(values = k_pal)
  
  return(p)
}