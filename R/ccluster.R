compute_slopes <- function(from, to, codebook, basis) {
  Beta <- do.call(cbind, codebook$coefs)
  x <- c(from, to)
  X <- basis(x)
  Y <- X %*% Beta
  slopes <- (Y[2, ] - Y[1, ]) / (to - from)
}

compute_feat1 <- function(codebook, basis) {
  compute_slopes(0, 5, codebook, basis)
}

compute_feat2 <- function(codebook, basis) {
  compute_slopes(5, 10, codebook, basis)
}

compute_feat3 <- function(codebook, basis) {
  compute_slopes(0, 15, codebook, basis)
}

compute_penalty <- function(curves, clusters, features, window) {
  w <- window / 2
  ncurves <- length(curves)
  penalties <- lapply(seq(ncurves), function(cx) {
    crv <- curves[[cx]]
    crv_feats <- crv[["features"]]    
    clust <- clusters[cx]
    clust_feats <- features[, clust]
    p <- abs(crv_feats - clust_feats) - w
    p[p < 0] <- 0
    sum(na.omit(p^2))
  })
  sum(do.call(c, penalties))
}

fksplines <- function(ksp, window) {
  algo <- ksp$algo
  basis <- ksp$basis
  cb <- ksp$codebook
  clusters <- ksp$clusters
  curves <- ksp$curves
  logl <- ksp$final_logl
  k <- ksp$k

  ## Split move---look for cluster with best dimension-wise improvement.
  for (cx in seq(k)) {
    
  }
}

#' Refine ksplines fit using feature constraints
#'
#' @export
cksplines <- function(ksp, window, temperatures, mc_len) {
  algo <- ksp$algo
  basis <- ksp$basis
  cb <- ksp$codebook
  clusters <- ksp$clusters
  curves <- ksp$curves
  logl <- ksp$final_logl
  k <- ksp$k
  features <- rbind(compute_feat1(cb, basis),
                    compute_feat2(cb, basis),
                    compute_feat3(cb, basis))
  penalty <- compute_penalty(curves, clusters, features, window)
  energy <- -logl + penalty

  ## For each temperature
  for (T in temperatures) {
    message(sprintf("Temp=%f", T))
    message(sprintf("LL=%f", logl))
    message(sprintf("Penalty=%f", penalty))
    message(sprintf("Energy=%f", energy))
    naccepts <- 0
    for (ix in seq(mc_len)) {
      prop_clusters <- sample(k, length(curves), replace=TRUE)
      prop_cb <- fit_codebook(k, curves, prop_clusters, algo)
      prop_logl <- full_logl(curves, prop_cb)
      prop_penalty <- compute_penalty(curves, prop_clusters, features, window)
      prop_energy <- -prop_logl + prop_penalty
      ## message(sprintf("Proposed=(%f)+(%f)=%f", -prop_logl, prop_penalty, prop_energy))
      if (prop_energy < energy) {
        clusters <- prop_clusters
        cb <- prop_cb
        logl <- prop_logl
        penalty <- prop_penalty
        energy <- prop_energy
        naccepts <- naccepts + 1
      } else {
        d <- prop_energy - energy
        p <- exp(- d / T)
        if (runif(1) < p) {
          clusters <- prop_clusters
          cb <- prop_cb
          logl <- prop_logl
          penalty <- prop_penalty
          energy <- prop_energy
          naccepts <- naccepts + 1
        }
      }
    }
    message(sprintf("AcceptPerc=%f", naccepts / mc_len))
  }

  ksp$cb <- cb
  ksp$clusters <- clusters
  ksp$final_logl <- logl
  ksp
}
