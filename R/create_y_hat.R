#' Multiply weighted adjacency matrix by control variables and resampled residuals
#'
#' @param beta Point estimates.
#' @param X Matrix of control variables.
#' @param par_dep_mat Weighted adjacency matrix.
#' @param sample_res Resampled residuals.
#' @return New dependent variable.
#' @keywords internal
#' @importFrom MASS ginv
#' @noRd
create_y_hat <- function(beta, X, par_dep_mat, sample_res){
  bX <- X %*% beta
  y_hat <- ginv(par_dep_mat) %*% (bX + sample_res)
  return(y_hat)
}
