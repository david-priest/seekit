# plotAbundancesA1.R — CATALYST-free rewrite promoted from dev/catalyst_quarantine/.
# Body verbatim from the project-vendored copy; only the CATALYST namespace
# shims (CATALYST:::.* internals, bare accessors, the asNamespace('CATALYST')
# hack) were rewritten to the package's .wl_* internals (R/wl_internals.R).
# 2026-06 seekit migration of the CMV CyTOF pipeline.
# plotAbundancesA1
# 
# Abundances plot variant A1 (lines 1021-1140 of source .Rmd).
# Migrated from CMV CyTOF Figures David.Rmd as part of repository reorganisation
# (see CMV_paper_analysis.Rmd / CMV_extra_analyses.Rmd / CMV_code_quarantine.Rmd).


plotAbundancesA1 <- function (x, 
                              k = "meta20", 
                              log = F, 
                              x_val = "day2", 
                              fill_by = NULL,
                              shape_by = NULL, 
                              title = " ",             # What meta to extract from sce colData
                              panel_spacing = 4,
                              miny = 0.01,
                              meta = c("sample_id","patient_id","condition"),
                              my_palette = NULL,
                              maxy = 100,
                              col_clust = TRUE,
                              k_pal = .wl_cluster_cols,
                              point_size = 1, 
                              lwidth = 0.4,
                              by = c("sample_id", "cluster_id"),
                              n_cols = 4, 
                              xaxis = "free_y", 
                              facet_ratio = 0.6,
                              excluded_clusters = NULL,
                              excluded_donors = NULL,
                              merging_col = NULL,
                              textsize = 12)
{
  
  
  
  # Use the merging column from colData if provided
  if (!is.null(merging_col)) {
    cluster_ids <- x[[k]]
  } else {
    k <- .wl_check_k(x, k)
    cluster_ids <- .wl_cluster_ids(x, k)
  }
  
  ns <- table(cluster_id = cluster_ids, sample_id = .wl_sample_ids(x))
  fq <- prop.table(ns, 2) * 100
  df <- as.data.frame(fq)
  m <- match(df$sample_id, x$sample_id)
  for (i in meta) df[[i]] <- x[[i]][m]
  
  dfout1 <<- df
  
  # Decide whether to plot on a log scale
  if (log == TRUE) {
    df$Freq <- df$Freq + 0.01
    ylabs <- "Proportion [%] (Log10 scale with 0.01 added to zero values)"
  }
  else {
    ylabs <- "Proportion [%]"
  }
  
  p <- ggplot(na.omit(df), aes_string(x = x_val, y = "Freq", fill = fill_by)) + 
    labs(x = NULL, y = ylabs) + 
    theme_bw() + 
    theme(
      panel.grid = element_blank(),
      text = element_text(size = textsize), 
      strip.text = element_text(size = textsize),
      strip.background = element_rect(fill = NA,color = NA),
      legend.text = element_text(size = textsize),
      aspect.ratio = facet_ratio,
      panel.spacing = unit(panel_spacing, "lines"),
      legend.title = element_text(size = textsize),
      axis.text = element_text(color = "black", size = textsize),
      axis.title = element_text(size = textsize),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = textsize),  # change if you want to include the x labels.
      #axis.text.x = element_blank(),
      axis.text.y = element_text(margin = margin(r = 10)),
      legend.key.height = unit(0.8, "lines"),
      axis.line = element_line(color = "black", linewidth = lwidth),  # Adjust the size here for axis line
      axis.ticks = element_line(color = "black",linewidth = lwidth),  # Adjust the size here for ticks
      panel.border = element_rect(color = "black", size = lwidth),  # Adjust the size here for panel border
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + 
    #geom_boxplot(outlier.colour = NA) +
    geom_boxplot(color = "black", position = position_dodge(), size = lwidth, alpha = 0.8, outlier.color = NA) +
    #geom_quasirandom(aes_string(shape = shape_by, group = fill_by), fill = "grey84",dodge.width = 0.75,width = 0.08, size = point_size, stroke = 0.5, shape = 21,show.legend = TRUE) + # , shape = 21
    geom_quasirandom(aes_string(shape = shape_by, group = fill_by, fill = fill_by),dodge.width = 0.75,width = 0.08, size = point_size, stroke = 0.5, shape = 21,show.legend = TRUE) + # , shape = 21
    #geom_line(aes_string(group = "patient_id"), alpha = 0.4, linewidth = 0.2) +
    scale_color_manual(values = my_palette) +
    scale_fill_manual(values = my_palette) +
    #facet_wrap2(~cluster_id, ncol = n_cols,scales = xaxis, axes = "y", remove_labels = "x") +
    facet_wrap2(~cluster_id, ncol = n_cols, scales = xaxis, axes = "all") +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold")) + labs(x = NULL, y = ylabs)
  
  dfout2 <<- df 
  
  if (log == TRUE) {
    if (xaxis == "free") {
      p + scale_y_continuous(trans='log10') + annotation_logticks(base = 10, sides = "l", outside = TRUE, size = lwidth) + coord_cartesian(clip = "off") + scale_size_area(max_size = 15)
    }
    else
    {
      p + scale_y_continuous(trans='log10',limits = c(miny,maxy),breaks=c(0.01,0.1,1,10,100),labels=c(0.01,0.1,1,10,100)) + annotation_logticks(base = 10, sides = "l", outside = TRUE, size = lwidth) + coord_cartesian(clip = "off") + scale_size_area(max_size = 15)
    }
  }
  else { # Log = FALSE
    if (xaxis == "free"){
      p + coord_cartesian(clip = "off") + scale_size_area(max_size = 15)  # plot without log free y axis and showing x axis ticks
    }
    else
    {
      p + scale_y_continuous(limits = c(0,maxy)) + coord_cartesian(clip = "off") + scale_size_area(max_size = 15)  # plot without log
    }
  }
  
} # end of function

