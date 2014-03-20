#' Learn spline codebook for sparse curves
#'
#' @export
#' 
ksplines <- function(curves, k, df, order, lambda = 0.1, maxiter = 10, normalize = FALSE) {
    all_x <- lapply(curves, "[[", "x")
    all_x <- do.call(c, all_x)
    basis_func <- get_basis(all_x, df, order)

    Omega <- smoothness_penalty(attr(basis_func, "knots"),
                                attr(basis_func, "bknots"),
                                attr(basis_func, "order"))
    Omega <- lambda * Omega
    
    for (ix in seq_along(curves)) {
        x <- curves[[ix]][["x"]]
        X <- basis_func(x)
        curves[[ix]][["X"]] <- X
    }
    
    fit_curves <- function(curves) {
        X <- do.call(rbind, lapply(curves, "[[", "X"))
        y <- do.call(c, lapply(curves, "[[", "y"))
        solve(crossprod(X) + Omega, crossprod(X, y))
    }

    fit_curves_norm <- function(curves) {
        all_X <- lapply(curves, function(curve) {
            n <- length(curve[["x"]])
            W <- (1 / sqrt(n)) * diag(n)
            W %*% curve[["X"]]
        })
        X <- do.call(rbind, all_X)

        all_y <- lapply(curves, function(curve) {
            n <- length(curve[["x"]])
            W <- (1 / sqrt(n)) * diag(n)
            as.numeric(W %*% curve[["y"]])
        })
        y <- do.call(c, all_y)

        solve(crossprod(X) + Omega, crossprod(X, y))
    }

    fit_codebook <- function(curves, cluster) {
        codebook <- matrix(0, nrow = df, ncol = k)
        
        for (cx in seq(k)) {
            subcurves <- curves[cluster == cx]
            
            if (length(subcurves) == 0) {
                codebook[, cx] <- 0
                next
            }

            if (normalize) {
                codebook[, cx] <- fit_curves_norm(subcurves)
            } else {
                codebook[, cx] <- fit_curves(subcurves)
            }
        }

        codebook
    }

    cluster <- sample(k, length(curves), replace = TRUE)
    codebook <- fit_codebook(curves, cluster)

    for (itx in seq(maxiter)) {
        for (ix in seq_along(curves)) {
            
            X <- curves[[ix]][["X"]]
            y <- curves[[ix]][["y"]]

            Yhat <- X %*% codebook
            r2 <- colSums((Yhat - y)^2)
            cluster[ix] <- which.min(r2)
        }

        codebook <- fit_codebook(curves, cluster)
    }

    new_ksplines(
        curves,
        cluster,
        codebook,
        basis_func)
}

#' Create new ksplines object
#'
new_ksplines <- function(curves, cluster, codebook, basis_func) {
    obj <- structure(list(), class = "ksplines")
    obj$curves <- curves
    obj$cluster <- cluster
    obj$codebook <- codebook
    obj$basis_func <- basis_func
    obj
}

