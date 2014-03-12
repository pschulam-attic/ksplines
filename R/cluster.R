#' Learn spline codebook for sparse curves
#'
#' @export
#' 
ksplines <- function(curves, k, df, degree, intercept = TRUE, lambda = 0.1, maxiter = 10) {
    all_x <- lapply(curves, "[[", "x")
    all_x <- do.call(c, all_x)
    basis_func <- choose_basis(all_x, df, degree, intercept)
    
    for (ix in seq_along(curves)) {
        x <- curves[[ix]][["x"]]
        X <- basis_func(x)
        curves[[ix]][["X"]] <- X
    }
    
    fit_curves <- function(curves) {
        X <- do.call(rbind, lapply(curves, "[[", "X"))
        y <- do.call(c, lapply(curves, "[[", "y"))
        solve(crossprod(X) + lambda * diag(df), crossprod(X, y))
    }

    fit_codebook <- function(curves, cluster) {
        codebook <- matrix(0, nrow = df, ncol = k)
        
        for (cx in seq(k)) {
            subcurves <- curves[cluster == cx]
            if (length(subcurves) == 0) {
                codebook[, cx] <- 0
                next
            }
            codebook[, cx] <- fit_curves(subcurves)
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

    fit <- structure(list(), class = "ksplines")
    fit$curves <- curves
    fit$cluster <- cluster
    fit$codebook <- codebook
    fit$basis_func <- basis_func
    fit
}

#' Choose an appropriate basis function for data 'x'
#' 
choose_basis <- function(x, df, degree, intercept) {
    B <- bs(x, df = df, degree = degree, intercept = intercept)
    knots <- attr(B, "knots")
    Boundary.knots <- attr(B, "Boundary.knots")

    function(x) {
        bs(x, knots = knots, Boundary.knots = Boundary.knots, intercept = intercept)
    }
}
