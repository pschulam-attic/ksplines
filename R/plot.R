#' @export
plot.ksplines <- function(fit, include_means = TRUE, print_immediately = TRUE) {
    curves <- fit$curves    
    curve_ids <- vapply(curves, "[[", integer(1), "id")
    curve_clusters <- data.frame(id = curve_ids, cluster = fit$cluster)

    data <- ldply(curves, function(curve) {
        data.frame(id = curve$id,
                   x = curve$x,
                   y = curve$y)
    })

    data <- merge(data, curve_clusters, by = "id")

    means <- data.frame()
    codebook <- fit$codebook
    xlow <- min(data$x)
    xhigh <- max(data$x)
    x <- seq(xlow, xhigh, length = 100)
    X <- fit$basis_func(x)

    for (cx in unique(data$cluster)) {
        mean_curve <- data.frame(cluster = cx,
                                 x = x,
                                 y = X %*% codebook[, cx])
        means <- rbind(means, mean_curve)
    }

    p <- ggplot() + theme_bw()
    p <- p + geom_point(aes(x = x, y = y), alpha = 0.25, data = data)
    p <- p + geom_line(aes(x = x, y = y), size = 1.5, color = "firebrick1", data = means)
    p <- p + facet_wrap(~ cluster)
    if (print_immediately) print(p)
    invisible(p)
}
