logsumexp <- function(x) {
  m <- max(x)
  log(sum(exp(x - m))) + m
}
