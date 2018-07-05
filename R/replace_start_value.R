#' Replace starting values with fitted values
#'
#' @param fit First object in the list of outcomes returned by net_dep. Passed by the function \code{boot}.
#' @param start List of starting values. Passed by the function \code{boot}.
#' @return New list of starting values.
#' @keywords internal
#' @noRd
replace_start_value <- function(fit, start){
  tmp <- coef(fit)
  for(i in 1:length(start)){
    start[[i]] <- as.numeric(tmp[names(tmp) %in% names(start)[i]])
  }
  start
}
