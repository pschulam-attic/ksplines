#' Choose an appropriate basis function for data 'x'
#' 
get_basis <- function(x, df, order) {
    B <- bs(x, df = df, degree = order - 1, intercept = TRUE)
    knots <- attr(B, "knots")
    Boundary.knots <- attr(B, "Boundary.knots")

    basis_func <- function(x) {
        bs(x, knots = knots, Boundary.knots = Boundary.knots, intercept = TRUE)
    }

    attr(basis_func, "knots") <- knots
    attr(basis_func, "bknots") <- Boundary.knots
    attr(basis_func, "order") <- order

    basis_func
}

#' Get smoothness penalty
#'
smoothness_penalty <- function(knots, bknots, order) {
    knots <- c(rep(bknots[1], order), knots, rep(bknots[2], order))
    sb <- SplineBasis(knots)
    OuterProdSecondDerivative(sb) / diff(bknots)
}
