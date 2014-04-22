logsumexp <- function(x) {
  m <- max(x)
  log(sum(exp(x - m))) + m
}

#' @export
estimate_features <- function(curve) {
  x <- curve[["x"]]
  y <- curve[["y"]]
  features <- rep(NA, 3)

  idx <- 0 <= x & x <= 5
  if (sum(idx) > 2) {
    x1 <- x[idx]
    y1 <- y[idx]
    lmfit <- lm(y1 ~ x1)
    features[1] <- coef(lmfit)[2]
  }

  idx <- 5 <= x & x <= 10
  if (sum(idx) > 2) {
    x2 <- x[idx]
    y2 <- y[idx]
    lmfit <- lm(y2 ~ x2)
    features[2] <- coef(lmfit)[2]
  }

  idx <- 0 <= x & x <= 15
  if (sum(idx) > 2) {
    x3 <- x[idx]
    y3 <- y[idx]
    lmfit <- lm(y3 ~ x3)
    features[3] <- coef(lmfit)[2]
  }

  curve[["features"]] <- features
  curve
}
