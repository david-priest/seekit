# plotClusterExprs2.R — CATALYST-free rewrite (.check_*/cluster_ids/.get_features/.agg -> .wl_*, hack removed).
# plotClusterExprs2
# 
# PlotExprHeatmap1 / cluster-expression plotting (formerly inline in CMV CyTOF Figures David.Rmd lines 32-100).
# Migrated from CMV CyTOF Figures David.Rmd as part of repository reorganisation
# (see CMV_paper_analysis.Rmd / CMV_extra_analyses.Rmd / CMV_code_quarantine.Rmd).

plotClusterExprs2 <- function (x, k = "meta20", features = "type", cluster_rows = TRUE, add_border = TRUE, text_size = 12, heatmap_palette = NULL, y_shift = 0.01, panel_spacing = 0.4, alpha_amount = 0.6) 
{
    .wl_check_sce(x, TRUE)
    k <- .wl_check_k(x, k)
    x$cluster_id <- .wl_cluster_ids(x, k)
    features <- .wl_get_features(x, features)
    
    ms <- t(.wl_agg(x[features, ], "cluster_id", "median"))
    
    if (cluster_rows) {
        d <- dist(ms, method = "euclidean")
        o <- hclust(d, method = "average")$order
    } else {
        o <- seq_along(levels(x$cluster_id))
    }
    
    cd <- colData(x)
    es <- assay(x[features, ], "exprs")
    df <- data.frame(t(es), cd, check.names = FALSE)
    df <- melt(df, id.vars = names(cd), variable.name = "antigen", 
               value.name = "expression")
    
    fq <- tabulate(x$cluster_id)/ncol(x)
    fq <- round(fq * 100, 2)
    names(fq) <- levels(x$cluster_id)
    cluster_levels <- rev(levels(x$cluster_id)[o])
    cluster_labels <- rev(names(fq)[o])
    df$cluster_id <- factor(df$cluster_id, levels = cluster_levels, labels = cluster_labels)
    
    # If a heatmap_palette is provided, use it to fill the colors, otherwise use red
    if (!is.null(heatmap_palette)) {
        cluster_colors <- heatmap_palette[levels(df$cluster_id)]
        scale_fill <- scale_fill_manual(values = cluster_colors)
    } else {
        scale_fill <- scale_fill_manual(values = rep("red", length(levels(df$cluster_id))))
    }
    
    # Define the base plot
    plot <- ggplot(df, aes(x = expression, y = cluster_id, fill = cluster_id)) +
        facet_wrap(~antigen, scales = "free_x", nrow = 1) + 
        geom_density_ridges(alpha = alpha_amount, rel_min_height = 0.01, scale = 1.1) + 
        scale_fill +
        scale_x_continuous(breaks = seq(0, 10, by = 2), limits = c(-0.5, NA), expand = c(0, 0)) +
        scale_y_discrete(expand = expand_scale(mult = c(y_shift, 0.2))) +
        theme(
              legend.position = "none", 
              strip.background = element_blank(), 
              strip.text = element_text(face = "bold", color = "black", size = text_size),
              panel.spacing = unit(panel_spacing, "lines"),
              plot.margin = unit(c(0, 1, 0, 0), "lines"),
              panel.background = element_rect(fill = "transparent", color = NA), # Make panel background transparent
              panel.border = if (add_border) element_rect(color = "grey43", fill = NA, size = 0.8) else element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(color = "black"), # Set x major gridlines as black
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              axis.line.x = element_blank(), # Hide main x-axis line
              axis.line.y = element_blank(), # Hide main y-axis line
              axis.text.x = element_text(color = "black"),
              axis.text.y = element_text(color = "black"),
              text = element_text(size = text_size))
    
    plot
}