#' Create network partitions with longitudinal data.
#'
#' @param x a vector of positive integer indexing the membership to one of the two parties.
#' @param G an object of class \code{Matrix} representing the social network.
#' @return a list of two objects: i) the matrix of connections within parties; ii) the matrix of connections between parties.
#' @keywords internal
#' @noRd
time_party_connection <- function(x, G, tt) {

  n <- unique(tt)

  partition_t <- lapply(n, function(y) x[tt == y])
  cong_mat <- lapply(partition_t, function(y) {
    tmp <- outer(y, y, "==")
    matrix(tmp * 1, length(y), length(y))
  } )
  within <- as.matrix(Reduce("bdiag", cong_mat))

  partition_t <- lapply(n, function(y) x[tt == y])
  cong_mat <- lapply(partition_t, function(y) {
    tmp <- outer(y, y, "==")
    tmp <- tmp - 1
    matrix(tmp * 1, length(y), length(y))
  } )
  between <- as.matrix(Reduce("bdiag", lapply(cong_mat, abs)))

  G_between <- G * between
  G_within <- G * within

  list(G_within, G_between)
}
