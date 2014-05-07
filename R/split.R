#' @export
num_feasible <- function(x, w, center=NULL) {
  x <- na.omit(x)
  if (is.null(center)) {
    center <- median(x)
  }
  feasible <- center - w / 2 <= x & x <= center + w / 2
  sum(feasible)
}

#' @import flexclust
#' @export
best_split <- function(curves, windows) {
  features <- do.call(rbind, lapply(curves, "[[", "features"))
  gain <- 0
  clusters <- NULL
  for (cx in seq(ncol(features))) {
    window <- windows[cx]
    f <- features[, cx]
    num_feats <- sum(!is.na(f))
    if (num_feats < 3) {
      next
    }
    nfeas <- num_feasible(f, window)
    x <- f
    x[is.na(x)] <- mean(x, na.rm=TRUE)
    kmfit <- kcca(x, 2, family = kccaFamily("kmedians"))
    f1 <- f[kmfit@cluster == 1]
    f2 <- f[kmfit@cluster == 2]
    split_nfeas <- num_feasible(f1, window) + num_feasible(f2, window)
    if (split_nfeas - nfeas > gain) {
      gain <- split_nfeas - nfeas
      clusters <- kmfit@cluster
    }
  }
  list(gain=gain, clusters=clusters)
}

#' @export
split_step <- function(ksp, window) {
  k <- ksp$k
  scores <- lapply(seq(k), function(cx) {
    subcurves <- ksp$curves[ksp$clusters == cx]
    best_split(subcurves, window)
  })
  gains <- vapply(scores, "[[", numeric(1), "gain")
  if (all(gains <= 0)) {
    message("No further splits")
    return(ksp)
  }
  best_split <- which.max(gains)
  best_clusters <- scores[[best_split]]$clusters
  if (is.null(best_clusters)) {
    message("No further splits.")
    return(ksp)
  }
  clusters <- numeric(length(ksp$clusters))
  clusters[ksp$clusters == best_split] <- best_clusters
  c1 <- which(clusters == 1)
  c2 <- which(clusters == 2)
  ksp <- ksplines_split(ksp, best_split, c1, c2)
  ksplines_sort(ksp)
}

#' @export
score_merge <- function(merge, ksp, windows) {
  cx1 <- merge[1]
  cx2 <- merge[2]
  sub1 <- ksp$curves[ksp$clusters == cx1]
  sub2 <- ksp$curves[ksp$clusters == cx2]
  feats1 <- do.call(rbind, lapply(sub1, "[[", "features"))
  feats2 <- do.call(rbind, lapply(sub2, "[[", "features"))
  
  nfeasible1 <- nfeasible2 <- 0
  for (cx in seq(along = windows)) {
    nfeasible1 <- nfeasible1 + num_feasible(feats1[, cx], windows[cx])
    nfeasible2 <- nfeasible2 + num_feasible(feats2[, cx], windows[cx])
  }
  
  feats <- rbind(feats1, feats2)
  nfeasible <- 0
  for (cx in seq(along = windows)) {
    nfeasible <- nfeasible + num_feasible(feats[, cx], windows[cx])
  }
  
  no_worse <- nfeasible - (nfeasible1 + nfeasible2) >= 0
  if (!no_worse) {
    list(logl_diff = -Inf, no_worse = no_worse)
  } else {
    ksp_merged <- ksplines_merge(ksp, cx1, cx2)
    logl_diff <- ksp_merged$final_logl - ksp$final_logl
    list(logl_diff = logl_diff, no_worse = no_worse)
  }
}

#' @export
merge_step <- function(ksp, windows) {
  k <- ksp$k
  merges <- combn(k, 2)
  merge_scores <- apply(merges, 2, score_merge, ksp, windows)
  no_worse <- vapply(merge_scores, "[[", logical(1), "no_worse")
  logl_diff <- vapply(merge_scores, "[[", numeric(1), "logl_diff")
  
  if (!any(no_worse & logl_diff > 0)) {
    ksp
  } else {
    merges <- merges[, no_worse]
    logl_diff <- logl_diff[no_worse]
    best_merge <- merges[, which.max(logl_diff)]
    ksplines_merge(ksp, best_merge[1], best_merge[2])
  }
}