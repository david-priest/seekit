# plot_colors
#
# Custom function extracted from CustomFunctions not Annotated New 250131.Rmd
# (lines 1658-1678). This is also exported by the CustomFunctionsSept25 R package;
# included here as a sourced script so the analysis .Rmd is self-contained.

plot_colors <- function(pal, title) {
  if (is.null(names(pal))) {
    names(pal) <- paste0(seq_along(pal), ": ", pal)
  }
  
  df <- data.frame(
    name = names(pal),
    color = pal,
    stringsAsFactors = FALSE
  )
  
  ggplot(df, aes(x = 1, y = factor(name, levels = rev(name)), fill = color)) +
    geom_tile(width = 2) +  # Increase the width of the tiles
    scale_fill_identity() +
    geom_text(aes(label = name), color = "black", size = 4, hjust = 0.5) +
    theme_void() +
    theme(legend.position = "none") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    coord_fixed(ratio = 1/10)  # Adjust the aspect ratio
}
