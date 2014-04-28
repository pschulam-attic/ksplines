codebook <- function(k) {
  cb <- structure(list(), class="codebook")
  cb$k <- k
  cb$probs <- rep(1, k) / k
  cb$coefs <- vector("list", k)
  cb$sigma <- 1
  cb
}

rearrange_codebook <- function(cb, p) {
  cb$probs <- cb$probs[p]
  cb$coefs <- cb$coefs[p]
  cb
}

#' @export
plot.codebook <- function(cb, basis, from, to) {
  x <- seq(from, to, length=100)
  X <- basis(x)
  p <- ggplot() + theme_bw()
  for (ix in seq(cb$k)) {
    y <- X %*% cb$coefs[[ix]]
    d <- data.frame(x=x, y=y)
    p <- p + geom_line(aes(x, y), data=d)
  }
  print(p)
  invisible(p)
}

codebook_posterior <- function(cb, X, y) {
  Coefs <- do.call(cbind, cb$coefs)
  Yhat <- X %*% Coefs
  res <- y - Yhat
  logl <- colSums(dnorm(res, sd=cb$sigma, log=TRUE))
  logjoint <- logl + log(cb$probs)
  logpost <- logjoint - logsumexp(logjoint)
  exp(logpost)
}

codebook_map <- function(cb, X, y) {
  posts <- codebook_posterior(cb, X, y)
  which.max(posts)
}

codebook_log_marginal <- function(cb, X, y) {
  Coefs <- do.call(cbind, cb$coefs)
  Yhat <- X %*% Coefs
  res <- y - Yhat
  logl <- colSums(dnorm(res, sd=cb$sigma, log=TRUE))
  logjoint <- logl + log(cb$probs)
  logsumexp(logjoint)
}
