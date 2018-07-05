#' Transform a matrix into an edge list
#' @param x an object of class \code{matrix}: adjacency matrix
#' @return an object of class \code{data.frame} edge list
#' @keywords internal
#' @noRd
mat_to_el <- function(x){
  delta <- row(x) - col(x)
  t(x)[which(delta != 0)]
}
