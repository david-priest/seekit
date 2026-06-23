# plotAbundanceDonut.R — CATALYST-free rewrite (.cluster_cols -> .wl_cluster_cols, namespace hack removed).
# plotAbundanceDonut
#
# Donut (polar-bar) plot of cluster abundances per group, with an optional
# nested "daughter" ring (e.g. a finer sub-clustering or isotype split drawn
# inside / on the edge of / across each master segment). Migrated faithfully
# from the iGCB nBass CyTOF analysis (241119 nBass Analysis 2) — drives the
# donut panels in paper Figs 1G, 2D and 2E.
#
# donut_type controls the daughter-ring geometry: "internal" (inside the master
# ring), "edges" (a separate outer ring), or "full" (full outer ring).
# average_across_samples = TRUE averages sample frequencies within each group
# and (via ei2) carries group-level metadata; FALSE plots one donut per sample.
#
# Cross-helper dep: ei2() (seekit). Runtime deps: dplyr, ggplot2, rlang,
# tidyr, ggnewscale.

plotAbundanceDonut <- function (x,
                                 k_abundances = "merging1",
                                 k_daughter = NULL,
                                 k_pal = .wl_cluster_cols,
                                 k_pal_daughter = NULL,
                                 group_by_order = NULL,
                                 sample_order = NULL,
                                 meta = c("sample_id", "patient_id", "condition"),
                                 group_by = "condition",
                                 shape_by = NULL,
                                 text_size = 16,
                                 n_cols = 4,
                                 rotang = 45,
                                 average_across_samples = TRUE,
                                 plot_ci = FALSE,
                                 cell_threshold = 0,
                                 drop_condition_if_any_sample_below_threshold = FALSE,
                                 plot_only_new_daughter = FALSE,
                                 donut_type = c("internal", "edges", "full")) {

  library(dplyr)
  library(ggplot2)
  library(rlang)
  library(tidyr)
  library(ggnewscale)

  donut_type <- match.arg(donut_type)

  #### MASTER DONUT DATA ####
  cluster_ids_abundances <- colData(x)[[k_abundances]]
  ns <- table(cluster_id = cluster_ids_abundances, sample_id = colData(x)$sample_id)
  dropped_samples <- character(0)
  dropped_conditions <- character(0)

  if (cell_threshold > 0) {
    sample_totals <- colSums(ns)
    if (drop_condition_if_any_sample_below_threshold && average_across_samples) {
      samp_df <- data.frame(sample_id = names(sample_totals), count = as.numeric(sample_totals),
                            stringsAsFactors = FALSE)
      meta_samp <- as.data.frame(colData(x)[, c("sample_id", group_by)], stringsAsFactors = FALSE)
      samp_df <- merge(samp_df, meta_samp, by = "sample_id")
      conditions_to_drop <- unique(samp_df[[group_by]][samp_df$count < cell_threshold])
      valid_conditions <- setdiff(unique(samp_df[[group_by]]), conditions_to_drop)
      valid_samples <- samp_df$sample_id[samp_df[[group_by]] %in% valid_conditions & samp_df$count >= cell_threshold]
      dropped_conditions <- conditions_to_drop
      dropped_samples <- setdiff(names(sample_totals), valid_samples)
    } else {
      valid_samples <- names(sample_totals)[sample_totals >= cell_threshold]
      dropped_samples <- names(sample_totals)[sample_totals < cell_threshold]
      meta_samp <- as.data.frame(colData(x)[, c("sample_id", group_by)], stringsAsFactors = FALSE)
      all_conditions <- unique(meta_samp[[group_by]])
      valid_cond <- unique(meta_samp[[group_by]][meta_samp$sample_id %in% valid_samples])
      dropped_conditions <- setdiff(all_conditions, valid_cond)
    }
    ns <- ns[, colnames(ns) %in% valid_samples, drop = FALSE]
    if (length(dropped_samples) > 0) {
      message("Dropped samples (total cell count < ", cell_threshold, "):")
      for (s in dropped_samples) message(" - ", s)
    }
    if (length(dropped_conditions) > 0) {
      message("Dropped conditions (no samples with total cell count >= ", cell_threshold, "):")
      for (c in dropped_conditions) message(" - ", c)
    }
  }

  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)
  m <- match(as.character(df$sample_id), as.character(colData(x)$sample_id))
  for (i in meta) {
    df[[i]] <- colData(x)[[i]][m]
  }

  if (average_across_samples) {
    df_avg <- df %>%
      group_by(cluster_id, !!sym(group_by)) %>%
      summarise(Freq = mean(Freq), .groups = "drop")
    df_ci <- df %>%
      group_by(cluster_id, !!sym(group_by)) %>%
      summarise(CI_lower = mean(Freq) - 1.96 * (sd(Freq) / sqrt(n())),
                CI_upper = mean(Freq) + 1.96 * (sd(Freq) / sqrt(n())), .groups = "drop")
    df2 <- left_join(df_avg, df_ci, by = c("cluster_id", group_by))
    df2 <- df2 %>%
      arrange(desc(cluster_id)) %>%
      group_by(!!sym(group_by)) %>%
      mutate(cumFreq = cumsum(Freq)) %>%
      ungroup()
    # [#6] BUGFIX: the original left_join'd ei2()'s per-sample table here, whose
    # `Freq` (count) column collided with df2's averaged `Freq` (and exploded
    # rows per sample) -> the default average_across_samples=TRUE path errored on
    # `mutate(Freq = Freq/sum(Freq))`. meta = group_by only re-adds the group
    # column df2 already has, so the join was a no-op; dropped it.
    df <- df2
    if (!is.null(group_by_order)) {
      df[[group_by]] <- factor(df[[group_by]], ordered = TRUE, levels = group_by_order)
    }
  } else {
    if (!is.null(sample_order)) {
      df$sample_id <- factor(df$sample_id, ordered = TRUE, levels = sample_order)
    }
  }

  # Re-normalize so the master percentages sum to exactly 100 WITHIN EACH FACET.
  # [#8] BUGFIX: the original always grouped by `group_by`, but with
  # average_across_samples = FALSE the facet is `sample_id` — grouping by the
  # condition then divided each sample's share by the number of samples in the
  # condition, so every donut only filled ~100/n_samples of the circle. Group by
  # the actual facet unit instead.
  norm_grp <- if (average_across_samples) group_by else "sample_id"
  df <- df %>%
    group_by(!!sym(norm_grp)) %>%
    mutate(Freq = Freq/sum(Freq)*100) %>%
    ungroup()

  #### DAUGHTER DONUT DATA ####
  daughter_data <- NULL
  if (!is.null(k_daughter) && !is.null(k_pal_daughter)) {
    # Retain each cell's master cluster (from k_abundances) and daughter value.
    daughter_df <- as.data.frame(colData(x)) %>%
      dplyr::select(sample_id, !!sym(group_by),
                    daughter = !!sym(k_daughter),
                    master = !!sym(k_abundances))
    daughter_counts <- daughter_df %>%
      group_by(sample_id, !!sym(group_by), master, daughter) %>%
      summarise(n = n(), .groups = "drop") %>%
      group_by(sample_id, !!sym(group_by)) %>%
      mutate(prop = n/sum(n)) %>%
      ungroup()
    if (average_across_samples) {
      daughter_data <- daughter_counts %>%
        group_by(!!sym(group_by), master, daughter) %>%
        summarise(Freq = mean(prop)*100, .groups = "drop")
    } else {
      daughter_data <- daughter_counts
      if (!is.null(sample_order))
        daughter_data$sample_id <- factor(daughter_data$sample_id, ordered = TRUE, levels = sample_order)
    }
    if (!is.null(group_by_order))
      daughter_data[[group_by]] <- factor(daughter_data[[group_by]], ordered = TRUE, levels = group_by_order)
    # Re-normalize so the daughter percentages sum to 100 WITHIN EACH FACET.
    # [#8] same facet-unit fix as the master ring above.
    dnorm_grp <- if (average_across_samples) group_by else "sample_id"
    daughter_data <- daughter_data %>%
      group_by(!!sym(dnorm_grp)) %>%
      mutate(Freq = Freq/sum(Freq)*100) %>%
      ungroup() %>%
      mutate(daughter = as.character(daughter),
             master = as.character(master))
  }

  #### PLOTTING ####
  # Set base x-coordinate and width based on donut_type
  if (donut_type == "internal") {
    master_x <- 4
    master_width <- 2
    daughter_x <- 4.5
    daughter_width <- 1
    xmin_new <- 4
    xmax_new <- 5
    x_limits <- c(2, 8)
  } else if (donut_type == "edges") {
    master_x <- 2
    master_width <- 1
    daughter_x <- 3
    daughter_width <- 1
    xmin_new <- 2.5
    xmax_new <- 3.2
    x_limits <- c(1, 4)
  } else { # full
    master_x <- 2
    master_width <- 1
    daughter_x <- 3.2
    daughter_width <- 1
    x_limits <- c(1, 4)
  }


  if (average_across_samples) {
    p <- ggplot(df, aes(x = master_x, y = Freq, fill = cluster_id)) +
      geom_bar(width = master_width, stat = "identity", color = "black", size = 0.2) +
      coord_polar(theta = "y", start = 0) +
      labs(x = NULL, y = "Proportion [%]") +
      theme_void() +
      theme(strip.text = element_text(face = "bold"),
            strip.background = element_rect(fill = NA, color = NA),
            text = element_text(size = text_size),
            legend.key.height = unit(0.8, "lines")) +
      facet_wrap(vars(!!sym(group_by)), ncol = n_cols)
  } else {
    p <- ggplot(df, aes(x = master_x, y = Freq, fill = cluster_id)) +
      geom_bar(width = master_width, stat = "identity", color = "black", size = 0.2) +
      coord_polar(theta = "y", start = 0) +
      labs(x = NULL, y = "Proportion [%]") +
      theme_void() +
      theme(strip.text = element_text(face = "bold"),
            strip.background = element_rect(fill = NA, color = NA),
            text = element_text(size = text_size),
            legend.key.height = unit(0.8, "lines")) +
      facet_wrap(vars(sample_id), ncol = n_cols)
  }
  p <- p + scale_fill_manual("Master", values = k_pal)

  if (!is.null(k_daughter) && !is.null(k_pal_daughter)) {
    p <- p + ggnewscale::new_scale_fill()
    if (!plot_only_new_daughter || donut_type == "full") {
      p <- p + geom_bar(data = daughter_data,
                        mapping = aes(x = daughter_x, y = Freq, fill = daughter),
                        stat = "identity", color = "black", size = 0.2, width = daughter_width,
                        inherit.aes = FALSE) +
        scale_fill_manual("Daughter", values = k_pal_daughter)
    } else {
      # Obtain master's boundaries from df.
      master_boundaries <- df %>%
        select(!!sym(group_by), master = cluster_id, master_Freq = Freq, cumFreq) %>%
        mutate(master_start = cumFreq - master_Freq)
      # Keep daughter rows where daughter label differs from the master.
      daughter_new <- daughter_data %>% filter(daughter != master)
      # Join the parent's start position.
      daughter_new <- left_join(daughter_new, master_boundaries, by = c(group_by, "master"))
      # For each (group, master), accumulate daughter segments so that the first starts exactly at master_start.
      daughter_new <- daughter_new %>%
        group_by_at(vars(!!sym(group_by), master)) %>%
        arrange(daughter) %>%
        mutate(cum_daughter = lag(cumsum(Freq), default = 0),
               d_ymin = master_start + cum_daughter,
               d_ymax = d_ymin + Freq) %>%
        ungroup()
      # Plot daughter segments with adjusted x-range
      p <- p + geom_rect(data = daughter_new,
                        mapping = aes(xmin = xmin_new, xmax = xmax_new,
                                      ymin = d_ymin, ymax = d_ymax, fill = daughter),
                        color = "black", size = 0.2, inherit.aes = FALSE) +
        scale_fill_manual("Daughter", values = k_pal_daughter)
    }
  }
    # Remove extra padding; adjust x-axis limits to accommodate larger donuts.
    # [#9] BUGFIX: limits = c(0,100) with the default oob = censor clips the top
    # stacked segment to NA when the cumulative floats just past 100 (e.g.
    # 100.0000001), leaving a wedge-gap at 12 o'clock. oob_squish pins it to 100
    # so the donut closes fully.
    p <- p + scale_x_continuous(limits = x_limits, expand = c(0,0)) +
           scale_y_continuous(limits = c(0,100), expand = c(0,0),
                              oob = scales::oob_squish)

  return(p)
}