fit <- function(algo, ...) {
  UseMethod("fit")
}

fit.leastsq <- function(lsq, X, y, w) {
  if (missing(w)) {
    w <- rep(1, length(y))
  }
  Xw <- X * w
  yw <- y * w
  A <- crossprod(Xw)
  b <- t(Xw) %*% yw
  as.numeric(solve(A, b))
}

fit.penleastsq <- function(plsq, X, y, w) {
  if (missing(w)) {
    w <- rep(1, length(y))
  }
  Xw <- X * w
  yw <- y * w
  A <- crossprod(Xw) + plsq$Omega
  b <- t(Xw) %*% yw
  as.numeric(solve(A, b))
}

leastsq <- function() {
  lsq <- structure(list(), class="leastsq")
  lsq
}

penleastsq <- function(Omega) {
  plsq <- structure(list(), class="penleastsq")
  plsq$Omega <- Omega
  plsq
}
