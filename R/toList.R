#' Function to transform a dataframe object in a list
#'
#' @param th an object of class \code{data.frame}
#' @return a list containing the vectors of the data frame \code{th}.
#' @keywords internal
#' @noRd
toList <- function(th) structure(as.list(th), names= names(th))
