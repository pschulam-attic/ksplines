#' @export
num_feasible <- function(x, w, center=NULL) {
  if (is.null(center)) {
    center <- mean(x)
  }
  feasible <- center - w / 2 <= x & x <= center + w / 2
  sum(feasible)
}

#' @export
split_scores <- function(curves, window) {
  features <- do.call(rbind, lapply(curves, "[[", "features"))
  gains <- numeric(ncol(features))
  for (cx in seq(ncol(features))) {
    f <- na.omit(features[, cx])
    if (length(f) < 3) {
      gains[cx] <- 0
      next
    }
    nfeas <- num_feasible(f, window)
    kmfit <- kmeans(f, 2)
    f1 <- f[kmfit$cluster == 1]
    f2 <- f[kmfit$cluster == 2]
    split_nfeas <- num_feasible(f1, window) + num_feasible(f2, window)
    gains[cx] <- split_nfeas - nfeas
  }
  gains
}

#' @export
split_step <- function(ksp, window) {
  k <- ksp$k
  scores <- lapply(seq(k), function(cx) {
    subcurves <- ksp$curves[ksp$clusters == cx]
    split_scores(subcurves, window)
  })
  scores <- do.call(rbind, scores)
  best_split <- arrayInd(which.max(scores), dim(scores))
  best_clust <- best_split[1, 1]
  best_feat <- best_split[1, 2]

  ksp$splits <- c(ksp$splits, best_clust)
  sub_idx <- ksp$clusters == best_clust
  subcurves <- ksp$curves[sub_idx]
  features <- vapply(subcurves, function(crv) crv$features[best_feat], numeric(1))
  features[is.na(features)] <- mean(features, na.rm=TRUE)
  kfit <- kmeans(features, 2)
  kclusters <- integer(length(ksp$clusters))
  kclusters[sub_idx] <- kfit$cluster
  ksp$clusters[sub_idx & kclusters == 2] <- ksp$k + 1
  ksp$k <- ksp$k + 1
  ksp$codebook$k <- ksp$k
  ksp$codebook$probs[ksp$k] <- ksp$codebook$probs[best_clust] / 2
  ksp$codebook$probs[best_clust] <- ksp$codebook$probs[ksp$k]

  sub1 <- ksp$curves[ksp$clusters == best_clust]
  ksp$codebook$coefs[[best_clust]] <- fit(ksp$algo, get_X(sub1), get_y(sub1))
  sub2 <- ksp$curves[ksp$clusters == ksp$k]
  ksp$codebook$coefs[[ksp$k]] <- fit(ksp$algo, get_X(sub2), get_y(sub2))
  
  return(ksp)
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
