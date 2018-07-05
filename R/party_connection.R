#' Create network partitions
#'
#' @param x a vector of positive integer indexing the membership to one of the two parties.
#' @param G an object of class \code{Matrix} representing the social network.
#' @return a list of two objects: i) the matrix of connections within parties; ii) the matrix of connections between parties.
#' @keywords internal
#' @noRd
party_connection <- function(x, G){

  tmp <- outer(x, x, "==")
  within <- matrix(tmp * 1, length(x), length(x))

  tmp <- tmp - 1
  between <- matrix(tmp * - 1, length(x), length(x))

  G_between <- G * between
  G_within <- G * within

  list(G_within, G_between)
}
