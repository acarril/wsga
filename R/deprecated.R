#' Deprecated function name
#'
#' `rddsga()` is the previous name of `wsga_rdd()`. It is retained as a
#' deprecated alias and will be removed in a future major release. New code
#' should call `wsga_rdd()` directly.
#'
#' @param ... Arguments forwarded to [wsga_rdd()].
#' @return The return value of [wsga_rdd()].
#' @seealso [wsga_rdd()]
#' @keywords internal
#' @export
rddsga <- function(...) {
  .Deprecated("wsga_rdd", package = "wsga")
  wsga_rdd(...)
}
