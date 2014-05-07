codebook <- function(k) {
  cb <- structure(list(), class="codebook")
  cb$k <- k
  cb$probs <- rep(1, k) / k
  cb$coefs <- vector("list", k)
  cb$sigma <- rep(1, k)
  cb
}

rearrange_codebook <- function(cb, p) {
  cb$probs <- cb$probs[p]
  cb$coefs <- cb$coefs[p]
  cb
}

codebook_reorder <- function(cb, perm) {
  cb$probs <- cb$probs[perm]
  cb$coefs <- cb$coefs[perm]
  cb
}

codebook_split <- function(cb, cx) {
  cb$k <- cb$k + 1
  cb$probs[cb$k] <- cb$probs[cx] / 2
  cb$probs[cx] <- cb$probs[cb$k]
  cb$coefs[[cb$k]] <- cb$coefs[[cx]]
  cb
}

codebook_merge <- function(cb, cx1, cx2) {
  cx <- min(cx1, cx2)
  cb$k <- cb$k - 1
  cb$probs[cx] <- sum(cb$probs[c(cx1, cx2)])
  cb$probs <- cb$probs[-max(cx1, cx2)]
  cb$coefs[[cx]] <- 0.5 * cb$coefs[[cx1]] + 0.5 * cb$coefs[[cx2]]
  cb$coefs <- cb$coefs[-max(cx1, cx2)]
  cb
}

codebook_multi_merge <- function(cb, cxs) {
  stopifnot(length(cxs) > 1)
  cx <- min(cxs)
  others <- cxs[-which.min(cxs)]
  cb$k <- cb$k - length(cxs) + 1
  cb$probs[cx] <- sum(cb$probs[cxs])
  cb$probs <- cb$probs[-others]

  v <- numeric(length(cb$coefs[[cx]]))
  for (ix in cxs) {
    v <- v + (1 / length(cxs)) * cb$coefs[[ix]]
  }
  cb$coefs[[cx]] <- v
  cb$coefs <- cb$coefs[-others]

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
  logl <- rowSums(dnorm(t(res), sd=cb$sigma, log=TRUE))
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
  logl <- rowSums(dnorm(t(res), sd=cb$sigma, log=TRUE))
  logjoint <- logl + log(cb$probs)
  logsumexp(logjoint)
}
