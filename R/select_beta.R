#' Select the names of coefficients associated to control variables in model specification
#'
#' @param specification Model specification
#' @return Vector of strings indicating coefficient names.
#' @keywords internal
#' @noRd
select_beta <- function(specification){

  beta <- all.vars(specification)
  sel_beta <- do.call("c", lapply(beta, function(x) grepl("beta_", x)))
  sel_alpha <- do.call("c", lapply(beta, function(x) grepl("alpha", x)))
  beta <- beta[which(sel_beta + sel_alpha == 1)]
  beta <- sub(" ", "", beta)
  beta <- sub(" ", "", beta)
  beta <- sub("\\(", "", beta)
  return(beta)
}

