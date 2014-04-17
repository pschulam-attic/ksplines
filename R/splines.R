#' @export
get_basis <- function(from, to, df, order) {
  if (df - order < 0) {
    stop("df must be >= order.")
  }
  num_iknots <- df - order
  total_knots = num_iknots + 2
  all_knots <- seq(from, to, length=total_knots)
  bknots <- all_knots[c(1, total_knots)]
  knots <- all_knots[-c(1, total_knots)]
  
  basis <- function(x) {
    bs(x, knots=knots, degree=order - 1,
       Boundary.knots=bknots, intercept=TRUE)
  }
  attr(basis, "knots") <- knots
  attr(basis, "bknots") <- bknots
  attr(basis, "order") <- order
  basis
}

#' @export
get_smoothness_penalty <- function(basis) {
  knots <- attr(basis, "knots")
  bknots <- attr(basis, "bknots")
  order <- attr(basis, "order")
  all_knots <- c(rep(bknots[1], order),
                 knots,
                 rep(bknots[2], order))
  
  sb <- SplineBasis(all_knots)
  OuterProdSecondDerivative(sb) / diff(bknots)
}
