#' @export
plot.ksplines <- function(fit) {
    curves <- fit$curves    
    curve_ids <- vapply(curves, "[[", integer(1), "id")
    curve_clusters <- data.frame(id = curve_ids, cluster = fit$cluster)

    data <- ldply(curves, function(curve) {
        data.frame(id = curve$id,
                   x = curve$x,
                   y = curve$y)
    })

    data <- merge(data, curve_clusters, by = "id")

    p <- ggplot(data) + theme_bw()
    p <- p + geom_point(aes(x = x, y = y), alpha = 0.25)
    p <- p + facet_wrap(~ cluster)
    p
}
