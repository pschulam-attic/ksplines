curve_logl <- function(crv, cb) {
  codebook_log_marginal(cb, crv[["X"]], crv[["y"]])
}

curve_map <- function(crv, cb) {
  codebook_map(cb, crv[["X"]], crv[["y"]])
}

full_logl <- function(curves, cb) {
  sum(vapply(curves, curve_logl, numeric(1), cb))
}

fit_codebook <- function(k, curves, clusters, algo) {
  cb <- codebook(k)
  res_sq <- 0
  npoints <- 0
  for (cx in seq(k)) {
    subcurves <- curves[clusters == cx]
    X <- do.call(rbind, lapply(subcurves, function(crv) {
      X <- crv[["X"]]
      X / sqrt(length(crv[["x"]]))
    }))
    y <- do.call(c, lapply(subcurves, function(crv) {
      y <- crv[["y"]]
      y / sqrt(length(crv[["x"]]))
    }))
    cb$probs[cx] <- length(subcurves)    
    cb$coefs[[cx]] <- fit(algo, X, y)
    res_sq <- res_sq + sum((y - X %*% cb$coefs[[cx]])^2)
    npoints <- npoints + length(y)
  }
  cb$probs <- cb$probs / length(curves)
  cb$sigma <- sqrt(res_sq / npoints)
  cb
}


#' Fit a spline mixture to curves
#' 
#' @export
ksplines <- function(curves, k, df, order, lambda, xrange=NULL) {
  if (is.null(xrange)) {
    all_x <- do.call(c, lapply(curves, "[[", "x"))
    xrange <- range(all_x)
  }
  
  basis <- get_basis(xrange[1], xrange[2], df, order)
  Omega <- lambda * get_smoothness_penalty(basis)
  algo <- penleastsq(Omega)
  for (ix in seq_along(curves)) {
    curves[[ix]][["X"]] <- basis(curves[[ix]][["x"]])
  }

  clusters <- sample(k, length(curves), replace=TRUE)
  cb <- fit_codebook(k, curves, clusters, algo)
  logl <- full_logl(curves, cb)
  repeat {
    old_logl <- logl
    clusters <- vapply(curves, curve_map, numeric(1), cb)
    cb <- fit_codebook(k, curves, clusters, algo)
    logl <- full_logl(curves, cb)
    if (logl - old_logl < 1e-2) {
      break
    }
  }

  ksp <- structure(list(), class="ksplines")
  ksp$k <- k
  ksp$codebook <- cb
  ksp$clusters <- clusters
  ksp$basis <- basis
  ksp$curves <- curves
  ksp
}
