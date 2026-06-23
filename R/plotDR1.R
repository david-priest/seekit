# plotDR1 is retired in favour of plotDR2 — plotDR2 is a superset (it absorbed
# plotDR1's unique args: border_width, plot_order, highlight_cluster). Kept as a
# silent forwarder so existing notebooks that call plotDR1 keep working unchanged.
# (Function wrapper, not `plotDR1 <- plotDR2`, so it's load-order-safe in the pkg.)
plotDR1 <- function(...) plotDR2(...)
