#' Create dyadic fixed effects
#'
#' @param x an object of class \code{matrix}: adjacency matrix.
#' @return an object of class \code{data.frame} containing fixed effect for both \emph{i} and \emph{j}.
#' @keywords internal
#' @noRd
create_dyadic_fe <- function(x){
  x <- as.matrix(x)
  ridx <- c(apply(x, 1, function(z) rep(z, length(x))))
  cidx <- rep(x, length(x))
  fixed <- cbind(ridx, cidx)
  fixed <- fixed[ - which(fixed[, 1] == fixed[, 2]), ]
  return(fixed)
}
