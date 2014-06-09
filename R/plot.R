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

#' @export
plot_curves_and_means <- function(ksp) {
  p <- ggplot() + theme_bw()
  p <- p + curve_layer(ksp, geom_point, alpha = 0.5)
  p <- p + mean_curve_layer(ksp, color = "red", size = 1.5)
  p + facet_wrap(~ class)
}

#' @export
plot_means_and_boot <- function(ksp) {
  p <- ggplot() + theme_bw()
  p <- p + bootstrap_layer(ksp, B = 100, alpha = 0.5)
  p <- p + mean_curve_layer(ksp, color = "red", size = 1.5)
  p + facet_wrap(~ class)
}

#' @export
curve_layer <- function(ksp, geom, ...) {
  data <- ldply(ksp$curves, function(crv) {
    data.frame(ptid = crv$id, x = crv$x, y = crv$y)
  })
  
  curve_ids <- vapply(ksp$curves, "[[", integer(1), "id")
  curve_classes <- data.frame(ptid = curve_ids, class = ksp$clusters)
  data <- merge(data, curve_classes, by = "ptid")
  
  geom(aes(x = x, y = y, group = ptid), data = data, ...)
}

#' @export
mean_curve_layer <- function(ksp, ...) {
  mean_curves <- data.frame()
  bknots <- attr(ksp$basis, "bknots")
  x_grid <- seq(bknots[1], bknots[2], length = 100)
  X_grid <- ksp$basis(x_grid)
  
  for (cx in seq(ksp$k)) {
    coefficients <- ksp$codebook$coefs[[cx]]
    yhat <- X_grid %*% coefficients
    mean_curve <- data.frame(x = x_grid, y = yhat, class = cx)
    mean_curves <- rbind(mean_curves, mean_curve)
  }
  
  geom_line(aes(x = x, y = y, group = class), data = mean_curves, ...)
}

#' @export
bootstrap_layer <- function(ksp, B = 100, ...) {
  algo <- ksp$algo
  k <- ksp$k
  
  bknots <- attr(ksp$basis, "bknots")
  x_grid <- seq(bknots[1], bknots[2], length = 100)
  X_grid <- ksp$basis(x_grid)
  
  ribbons <- data.frame()
  
  for (cx in seq(k)) {
    subcurves <- ksp$curves[ksp$clusters == cx]
    boot_coef <- bootstrap_coefficients(subcurves, algo, B)
    boot_means <- X_grid %*% boot_coef
    lower_curve <- apply(boot_means, 1, function(x) quantile(x, 0.025))
    upper_curve <- apply(boot_means, 1, function(x) quantile(x, 0.975))
    curves <- data.frame(x = x_grid, lower = lower_curve, upper = upper_curve, class = cx)
    ribbons <- rbind(ribbons, curves)
  }
  
  geom_ribbon(aes(x = x, ymin = lower, ymax = upper), data = ribbons, ...)
}