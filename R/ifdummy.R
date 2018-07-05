#' Create a dyadic dummy variable using continous variables#'
#' @param x a numeric vector: characteristics of agent \emph{i}.
#' @param y a numeric vector: characteristics of agent \emph{j}.
#' @return a dyadic dummy variable which takes 1 if both agents have the same characteristic, and 0 otherwise
#' @keywords internal
#' @noRd
ifdummy <- function(x, y){
  x <- as.numeric(as.character(x))
  y <- as.numeric(as.character(y))
  res <- x - y
  ifelse(res == 0, 1, 0)
}
