#' Create network partitions
#'
#' @param x a vector of positive integer indexing the membership to one of the two parties.
#' @param G an object of class \code{Matrix} representing the social network.
#' @return a list of two objects: i) the matrix of connections within parties; ii) the matrix of connections between parties.
#' @keywords internal
#' @noRd
party_connection_split <- function(x, G){

  res <- party_connection(x, G)
  G_within <- res[[1]]
  G_between <- res[[2]]

  G_within_0 <- G_within
  G_within_0[which(x %in% max(x)), ] <- 0
  G_within_0[, which(x %in% max(x))] <- 0

  G_within_1 <- G_within
  G_within_1[which(x %in% min(x)), ] <- 0
  G_within_1[, which(x %in% min(x))] <- 0

  G_between_01 <- G_between
  G_between_01[, which(x %in% min(x))] <- 0

  G_between_10 <- G_between
  G_between_10[, which(x %in% max(x))] <- 0

  list(G_within_0, G_within_1, G_between_01, G_between_10)
}
