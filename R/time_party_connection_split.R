#' Create network partitions with longitudinal data.
#'
#' @param x a vector of positive integer indexing the membership to one of the two parties.
#' @param G an object of class \code{Matrix} representing the social network.
#' @return a list of two objects: i) the matrix of connections within parties; ii) the matrix of connections between parties.
#' @keywords internal
#' @noRd
time_party_connection_split <- function(x, G, tt) {

  res <- time_party_connection(x, G, tt)
  G_within <- res[[1]]
  G_between <- res[[2]]

  G_within_0 <- G_within
  G_within_0[x == max(x), ] <- 0
  G_within_0[, x == max(x)] <- 0

  G_within_1 <- G_within
  G_within_1[x == min(x), ] <- 0
  G_within_1[, x == min(x)] <- 0

  G_between_01 <- G_between
  G_between_01[, min(x)] <- 0

  G_between_10 <- G_between
  G_between_10[, max(x)] <- 0

  list(G_within_0, G_within_1, G_between_01, G_between_10)
}
