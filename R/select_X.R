#' Select the names of control variables in model specification
#'
#' @param specification Model specification
#' @return Vector of strings indicating variable names.
#' @keywords internal
#' @noRd
select_X <- function(specification){
  X <- as.character(specification[[3]])
  X <- as.character(X[[3]])
  X <- strsplit(X, "\\+")
  X <- as.character(sapply(X[[1]], function(x) strsplit(x, "\\*")[[1]])[2, ])
  X <- sub(" ", "", X)
  X <- sub(" ", "", X)
  X <- sub("\\)", "", X)
  return(X)
}
