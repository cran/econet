#' Extract the adjacency matrices from the nls or mle2 object.
#'
#' @param fit an object of class \code{nls} or \code{mle2}
#' @param type string. "lim" refers to the linear-in-means model. "het" implements the model with heterogenous peer effects in Battaglini, Patacchini (2018). "het_l" and "het_r" uses the models for heterogenous peer effects in Battaglini, Leone Sciabolazza, Patacchini (2018). "par" implements the model with network partitions by Battaglini, Leone Sciabolazza, Patacchini (2018).
#' @return an object of class \code{matrix}: adjacency matrix.
#' @keywords internal
#' @noRd
dep_mat <- function(fit, type){
  switch(type,
         lim = fit$data[["G"]],
         het = list(fit$data[["G"]], fit$data[["G_heterogeneity"]]),
         het_l = list(fit$data[["G"]], fit$data[["G_heterogeneity"]]),
         het_r = list(fit$data[["G"]], fit$data[["G_heterogeneity"]]),
         par = list(fit$data[["G_within"]], fit$data[["G_between"]])
  )
}
