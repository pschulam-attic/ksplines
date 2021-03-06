#' @export
feature_ksplines <- function(curves, k, df, order, feature_fn,
                             lambda = 0.1,
                             maxiter = 10,
                             normalize = FALSE) {
    
    all_x <- lapply(curves, "[[", "x")
    all_x <- do.call(c, all_x)
    basis_func <- get_basis(all_x, df, order)

    Omega <- smoothness_penalty(attr(basis_func, "knots"),
                                attr(basis_func, "bknots"),
                                attr(basis_func, "order"))
    Omega <- lambda  * Omega

    for (ix in seq_along(curves)) {
        x <- curves[[ix]][["x"]]
        X <- basis_func(x)
        curves[[ix]][["X"]] <- X
    }

    solve_lsq <- function(X, y) {
        obj_fn <- function(coefs) {
            r2 <- (y - X %*% coefs)^2
            penalty <- as.numeric(coefs %*% Omega %*% coefs)
            sum(r2) + penalty
        }

        grad_fn <- function(coefs) {
            r <- as.numeric(y - X %*% coefs)
            penalty <- 2 * Omega %*% coefs
            gr <- -2 * colSums(X * r) + as.numeric(penalty)
            gr <- gr / sqrt(sum(gr^2))
            gr
        }

        betas <- optim(rep(0, df), obj_fn, grad_fn, method = "BFGS")$par
        betas
    }

    fit_curves <- function(curves) {
        X <- do.call(rbind, lapply(curves, "[[", "X"))
        y <- do.call(c, lapply(curves, "[[", "y"))
        solve_lsq(X, y)
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

        solve_lsq(X, y)
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

### TODOS
# 1. Update objective function to reflect feature funcs.
# 2. Update gradient function to reflect feature funcs.
# 3. Update class reassignment to reflect feature funcs.
# 
