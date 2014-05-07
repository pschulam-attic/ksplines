#' @export
plot.ksplines <- function(ksp) {
  p <- plot_ksplines(ksp)
  print(p)
}

#' @export
plot_ksplines <- function(ksp, facet=TRUE) {
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
  means <- do.call(rbind, means[ksp$codebook$probs > 0])

  p <- ggplot() + theme_bw()
  p <- p + geom_line(aes(x=x, y=y, group=id), alpha=0.25, data=data)
  p <- p + geom_line(aes(x=x, y=y), alpha = 0.25, size=1, color="firebrick1", data=means)
  if (facet) {
    p <- p + facet_wrap(~ cluster)
  }
  p
}

#' @export
plot_curves <- function(curves) {
  data <- ldply(curves, function(crv) {
    data.frame(ptid=crv$id, x=crv$x, y=crv$y)
  })
  p <- ggplot(data) + theme_bw()
  p <- p + geom_line(aes(x=x, y=y, group=ptid), alpha=0.5)
  p
}
