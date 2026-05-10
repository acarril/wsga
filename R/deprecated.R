#' Deprecated function name
#'
#' `rddsga()` is the previous name of `wsga()`. It is retained as a deprecated
#' alias and will be removed in a future major release. New code should call
#' `wsga()` directly.
#'
#' @param ... Arguments forwarded to [wsga()].
#' @return The return value of [wsga()].
#' @seealso [wsga()]
#' @keywords internal
#' @export
rddsga <- function(...) {
  .Deprecated("wsga", package = "wsga")
  wsga(...)
}
