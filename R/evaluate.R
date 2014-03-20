#' Sample observed curves from a cluster
#'
#' @export
#'
sample_curves <- function(ksplines, cluster, n) {
    curves <- ksplines$curves
    curves <- curves[ksplines$cluster == cluster]
    samp <- sample(curves, min(length(curves), n))
    
    samp <- lapply(samp, function(curve) {
        df <- data.frame(id = curve$id, x = curve$x, y = curve$y)
    })
    
    samp <- do.call(rbind, samp)
    rownames(samp) <- NULL
    samp
}

sample_curves_plot <- function(ksplines, cluster, n) {
    samp <- sample_curves(ksplines, cluster, n)
    p <- ggplot(samp) + theme_bw()
    p <- p + geom_line(aes(x, y, group = id))
    p
}
