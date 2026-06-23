# plotAbundanceScatter — migrated from CyTOF nBass_helpers.R into seekit (CATALYST-free).
# 2026-06-10: lifted verbatim, de-CATALYST'd (.wl_* internals, namespace hack removed).

plotAbundanceScatter <- function(x, k_abundances = "merging1",
                                 meta = c("sample_id", "patient_id", "condition", "exp"),
                                 facet_by = NULL, text_size = 16,
                                 k_pal = .wl_cluster_cols, show_r_squared = TRUE) {
  cluster_ids_abundances <- colData(x)[[k_abundances]]
  ns <- table(cluster_id = cluster_ids_abundances, sample_id = colData(x)$sample_id)
  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)
  m <- match(df$sample_id, colData(x)$sample_id)
  for (i in meta) df[[i]] <- colData(x)[[i]][m]

  replicates <- unique(df$exp)
  if (length(replicates) != 2) {
    stop("The data must contain exactly two replicates in the 'exp' column.")
  }
  df <- df %>% dplyr::filter(exp %in% replicates)
  df1 <- df %>% dplyr::filter(exp == replicates[1]) %>% dplyr::rename(Freq1 = Freq)
  df2 <- df %>% dplyr::filter(exp == replicates[2]) %>% dplyr::rename(Freq2 = Freq)
  df_wide <- dplyr::left_join(df1, df2, by = c("patient_id", "cluster_id", "day"))

  p <- ggplot(df_wide, aes(x = Freq1, y = Freq2, color = cluster_id)) +
    geom_point(aes(shape = day), size = 3) +
    labs(x = paste("Replicate", replicates[1]), y = paste("Replicate", replicates[2]), color = "Cluster ID") +
    theme_bw() +
    theme(panel.grid = element_blank(), strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey88", color = NA), aspect.ratio = 1,
          axis.text = element_text(color = "black"), text = element_text(size = text_size),
          panel.spacing = unit(1, "lines"), legend.key.height = unit(0.8, "lines")) +
    scale_color_manual(values = k_pal) + xlim(0, 100) + ylim(0, 100)

  if (!is.null(facet_by)) p <- p + facet_wrap(vars(!!sym(facet_by)), scales = "free_x")

  if (show_r_squared && requireNamespace("ggpubr", quietly = TRUE)) {
    p <- p + ggpubr::stat_cor(aes(x = Freq1, y = Freq2), method = "spearman",
                              label.x.npc = "left", label.y.npc = "top", size = text_size / 4,
                              cor.coef.name = "rho", cor.method = "spearman",
                              p.accuracy = 0.001, label.sep = "\n", inherit.aes = FALSE)
  }
  return(p)
}
