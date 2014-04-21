#' Fit a featurized spline mixture to curves
#'
#' @export
cksplines <- function(curves, k, df, order, lambda, xrange=NULL) {
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

  ksp <- ksplines(curves, k, df, order, lambda, xrange)
  clusters <- ksp$clusters
  cb <- ksp$codebook
  logl <- ksp$final_logl
  

  ## For each temperature

  ## Make some number of moves
}
