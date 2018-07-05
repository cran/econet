#' Function to resample residuals
#'
#' @param res Residuals of the model.
#' @param par_dep_mat Weighted adjacency matrix.
#' @param group \code{NULL} or vector of positive integers specifying the indices for resampling within groups.
#' @return Vector of bootstrapped residuals.
#' @keywords internal
#' @noRd
sample_residuals <- function(res, par_dep_mat, group){

  r_mat <- par_dep_mat %*% res

  if (is.null(group)) {
    r_mat <- sample(r_mat, replace = T)
  } else {
    id_group <- unique(group)
    for (i in 1:length(id_group)) {
      sel <- which(group == id_group[i])
      r_mat[sel] <- ifelse(length(sel) == 1, r_mat[sel], sample(r_mat[sel],
                                                                replace = T))
    }
  }

  return(r_mat)

}
