#' Create weighted adjacency matrix
#'
#' @param par_dep Value of the estimated decay-parameter.
#' @param G Adjacency Matrix.
#' @param I Identity Matrix.
#' @param type hypothesis: e.g. \code{"lim"}, \code{"het"}, \code{"het_l"}, \code{"het_r"}, \code{"par"}.
#' @return Weighted adjacency matrix.
#' @keywords internal
#' @noRd
par_dep_mat <- function(par_dep, G, I, type){
  switch(type,
         lim = (I - par_dep * G),
         het = (I - G[[1]] %*% (par_dep[1] * I + par_dep[2] * G[[2]])),
         het_l = (I - G[[1]] %*% (par_dep[1] * I - par_dep[2] * G[[2]])),
         het_r = (I - G[[1]] %*% (par_dep[1] * I - par_dep[2] * G[[2]])),
         par = (I - par_dep[1] * G[[1]] - par_dep[2] * G[[2]]),
         par_split_btw = (I - par_dep[1] * G[[1]] - par_dep[2] * G[[2]] - par_dep[3] * G[[3]]),
         par_split_with = (I - par_dep[1] * G[[1]] - par_dep[2] * G[[2]] - par_dep[3] * G[[3]]),
         par_split_with_btw = (I - par_dep[1] * G[[1]] - par_dep[2] * G[[2]] - par_dep[3] * G[[3]] - par_dep[4] * G[[4]])
  )
}
