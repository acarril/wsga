#' Compute kernel weights
#'
#' @param x Numeric vector: running variable (centered at cutoff = 0 expected).
#' @param bwidth Positive numeric: half-bandwidth.
#' @param kernel Character: `"uniform"` (default), `"triangular"`, or `"epanechnikov"`.
#' @param obs_weights Optional numeric vector of pre-existing observation weights.
#'   Multiplied into the kernel weights element-wise.
#' @return Numeric vector of kernel weights (0 outside bandwidth).
#' @noRd
kernel_weights <- function(x, bwidth, kernel = c("uniform", "triangular", "epanechnikov"),
                           obs_weights = NULL) {
  kernel <- match.arg(kernel)
  ax <- abs(x)

  w <- switch(kernel,
    uniform      = as.numeric(ax < bwidth),
    triangular   = pmax(0, bwidth - ax),
    epanechnikov = pmax(0, 3 / 4 * (bwidth^2 - ax^2))
  )

  if (!is.null(obs_weights)) w <- w * obs_weights
  w
}
