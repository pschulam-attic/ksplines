curve_logl <- function(crv, cb) {
  codebook_log_marginal(cb, crv[["X"]], crv[["y"]])
}

curve_map <- function(crv, cb) {
  codebook_map(cb, crv[["X"]], crv[["y"]])
}

full_logl <- function(curves, cb) {
  sum(vapply(curves, curve_logl, numeric(1), cb))
}

get_X <- function(curves) {
  do.call(rbind, lapply(curves, function(crv) {
    X <- crv[["X"]]
    X / sqrt(length(crv[["x"]]))
  }))
}

get_y <- function(curves) {
  do.call(c, lapply(curves, function(crv) {
    y <- crv[["y"]]
    y / sqrt(length(crv[["x"]]))
  }))
}

fit_sigma <- function(curves, clusters, codebook) {
  k <- codebook$k
  sigma <- numeric(k)
  for (cx in seq(k)) {
    subcurves <- curves[clusters == cx]
    X <- get_X(subcurves)
    y <- get_y(subcurves)
    yhat <- X %*% codebook$coefs[[cx]]
    res_sq <- sum((y - yhat)^2)
    sigma[cx] <- sqrt(res_sq / length(subcurves))
  }
  sigma
}

fit_codebook <- function(k, curves, clusters, algo, one_sd=FALSE) {
  df <- ncol(curves[[1]][["X"]])
  cb <- codebook(k)
  res_sq <- 0
  for (cx in seq(k)) {
    subcurves <- curves[clusters == cx]
    if (length(subcurves) == 0) {
      cb$probs[cx] <- 0
      cb$coefs[[cx]] <- numeric(df)
      next
    }
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

    if (one_sd) {
      res_sq <- res_sq + sum((y - X %*% cb$coefs[[cx]])^2)
    } else {
      res_sq <- sum((y - X %*% cb$coefs[[cx]])^2)
      cb$sigma[cx] <- sqrt(res_sq / length(subcurves))
    }
  }
  
  cb$probs <- cb$probs / length(curves)
  if (one_sd) {
    cb$sigma[1:k] <- sqrt(res_sq / length(curves))
  }
  cb
}


#' Fit a spline mixture to curves
#' 
#' @export
ksplines <- function(curves, k, df, order, lambda,
                     xrange=NULL, tol=1e-2, verbose=FALSE, seed = 1) {
  
  set.seed(1)
  
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
  cb <- fit_codebook(k, curves, clusters, algo, one_sd = TRUE)
  logl <- full_logl(curves, cb)
  if (verbose) message(sprintf("LL=%f", logl))
  repeat {
    old_logl <- logl
    clusters <- vapply(curves, curve_map, numeric(1), cb)
    cb <- fit_codebook(k, curves, clusters, algo)
    logl <- full_logl(curves, cb)
    if (verbose) message(sprintf("LL=%f", logl))
    if (logl - old_logl < tol) {
      break
    }
  }

  ksp <- structure(list(), class="ksplines")
  ksp$k <- k
  ksp$codebook <- cb
  ksp$clusters <- clusters
  ksp$basis <- basis
  ksp$algo <- algo
  ksp$curves <- curves
  ksp$final_logl <- logl
  ksp$merges <- list()
  ksplines_sort(ksp)
}

#' @export
sksplines <- function(nsplits, window, ...) {
  ksp <- ksplines(...)
  max_k <- ksp$k + nsplits
  while (ksp$k < max_k) {
    prev_k <- ksp$k
    ksp <- split_step(ksp, window)
    if (ksp$k == prev_k) {
      break
    }
  }
  ksp
}

#' @export
smksplines <- function(nsplits, window, ...) {
  ksp <- ksplines(...)
  max_k <- ksp$k + nsplits
  
  while (ksp$k < max_k) {
    prev_k <- ksp$k
    ksp <- split_step(ksp, window)
    message("Split finished...")
    
    if (ksp$k == prev_k) {
      repeat {
        old_k <- ksp$k
        ksp <- merge_step(ksp, window)
        if (ksp$k == old_k) {
          break
        }
      }
      message(sprintf("K=%d", ksp$k))
      break
    } else {
      repeat {
        old_k <- ksp$k
        ksp <- merge_step(ksp, window)
        if (ksp$k == old_k) {
          break
        }
      }
      message(sprintf("K=%d", ksp$k))
    }
  }
  ksp
}

#' @export
ksplines_split <- function(ksp, cx, c1, c2) {
  ksp$k <- ksp$k + 1
  ksp$codebook <- codebook_split(ksp$codebook, cx)
  ksp$clusters[c2] <- ksp$k

  X <- get_X(ksp$curves[c1])
  y <- get_y(ksp$curves[c1])
  ksp$codebook$coefs[[cx]] <- fit(ksp$algo, X, y)
  ksp$codebook$probs[cx] <- length(c1) / length(ksp$curves)  

  X <- get_X(ksp$curves[c2])
  y <- get_y(ksp$curves[c2])
  ksp$codebook$coefs[[ksp$k]] <- fit(ksp$algo, X, y)
  ksp$codebook$probs[ksp$k] <- length(c2) / length(ksp$curves)
  
  ksp$codebook$sigma <- fit_sigma(ksp$curves, ksp$clusters, ksp$codebook)
  ksp$final_logl <- full_logl(ksp$curves, ksp$codebook)

  ksp
}

#' @export
ksplines_merge <- function(ksp, cx1, cx2) {
  cx <- min(cx1, cx2)
  merge_move <- list(cx, which(ksp$clusters == cx1), which(ksp$clusters == cx2))
  ksp$k <- ksp$k - 1
  ksp$codebook <- codebook_merge(ksp$codebook, cx1, cx2)
  ksp$clusters[ksp$clusters == max(cx1, cx2)] <- cx
  others <- ksp$clusters > max(cx1, cx2)
  ksp$clusters[others] <- ksp$clusters[others] - 1

  X <- get_X(ksp$curves[ksp$clusters == cx])
  y <- get_y(ksp$curves[ksp$clusters == cx])
  ksp$codebook$coefs[[cx]] <- fit(ksp$algo, X, y)
  
  ksp$codebook$sigma <- fit_sigma(ksp$curves, ksp$clusters, ksp$codebook)
  ksp$final_logl <- full_logl(ksp$curves, ksp$codebook)

  nmerges <- length(ksp$merges)
  ksp$merges[[nmerges + 1]] <- merge_move
  ksp
}

#' @export
ksplines_multi_merge <- function(ksp, cxs) {
  cx <- min(cxs)
  others <- cxs[-which.min(cxs)]
  ksp$k <- ksp$k - length(cxs) + 1
  ksp$codebook <- codebook_multi_merge(ksp$codebook, cxs)
  ksp$clusters[ksp$clusters %in% others] <- cx

  for (ix in sort(others, decreasing=TRUE)) {
    idxs <- ksp$clusters > ix
    ksp$clusters[idxs] <- ksp$clusters[idxs] - 1
  }

  X <- get_X(ksp$curves[ksp$clusters == cx])
  y <- get_y(ksp$curves[ksp$clusters == cx])
  ksp$codebook$coefs[[cx]] <- fit(ksp$algo, X, y)

  ksp
}

#' @export
ksplines_unmerge <- function(ksp) {
  if (is.null(ksp$merges)) {
    return(ksp)
  }

  nmerges <- length(ksp$merges)
  last_merge <- ksp$merges[[nmerges]]
  ksp <- ksplines_split(ksp, last_merge[[1]], last_merge[[2]], last_merge[[3]])
  ksp$merges <- ksp$merges[-nmerges]
  ksplines_sort(ksp)
}

#' @export
ksplines_reorder <- function(ksp, perm) {
  ksp$codebook <- codebook_reorder(ksp$codebook, perm)
  clust_map <- seq(ksp$k)
  clust_map[perm] <- seq(ksp$k)
  ksp$clusters <- clust_map[ksp$clusters]
  ksp
}

#' @export
ksplines_sort <- function(ksp, decreasing=TRUE) {
  perm <- order(ksp$codebook$probs, decreasing=decreasing)
  ksplines_reorder(ksp, perm)
}
