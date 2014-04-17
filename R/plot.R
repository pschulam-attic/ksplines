#' @export
plot.ksplines <- function(ksp, include_means=TRUE) {
  curves <- ksp$curves
  ids <- vapply(curves, "[[", integer(1), "id")
  clusts <- data.frame(id=ids, cluster=ksp$clusters)
  data <- ldply(curves, function(crv) data.frame(id=crv$id, x=crv$x, y=crv$y))
  data <- merge(data, clusts, by="id")

  from <- min(data$x)
  to <- max(data$x)
  x <- seq(from, to, length=100)
  X <- ksp$basis(x)
  means <- lapply(seq(ksp$k), function(cx) {
    y <- X %*% ksp$codebook$coefs[[cx]]
    data.frame(x=x, y=y, cluster=cx)
  })
  means <- do.call(rbind, means)

  p <- ggplot() + theme_bw()
  p <- p + geom_point(aes(x=x, y=y), alpha=0.25, data=data)
  p <- p + geom_line(aes(x=x, y=y), size=1.5, color="firebrick1", data=means)
  p <- p + facet_wrap(~ cluster)
  print(p)
}

