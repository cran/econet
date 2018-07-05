#' Compute the Weighted Katz-Bonacich centrality
#'
#' @param second_step an object of class "nls" or "mle2"
#' @param hypothesis string. "lim" refers to the linear-in-means model. "het" implements the model with heterogenous peer effects in Battaglini, Patacchini (2018). "het_l" and "het_r" uses the models for heterogenous peer effects in Battaglini, Leone Sciabolazza, Patacchini (2018). "par" implements the model with network partitions by Battaglini, Leone Sciabolazza, Patacchini (2018).
#' @param I an object of class \code{matrix}. Identity matrix.
#' @param G an object of class \code{matrix} representing the social network.
#' @param e a list of objects containing the data used in the regression.
#' @return a numeric vector: Weighted Katz-Bonacich centrality.
#' @keywords internal
#' @importFrom stats coef
#' @noRd
parameter_dependent_centrality <- function(second_step, hypothesis, I, G, e) {
  if (hypothesis == "lim") {
    phi <- coef(second_step)[names(coef(second_step)) %in% "phi"]
    centrality <- solve_block(I - phi * G) %*% e[["Ones"]]
  } else if (hypothesis == "het") {
    phi <- coef(second_step)[names(coef(second_step)) %in% "phi"]
    gamma <- coef(second_step)[names(coef(second_step)) %in% "phi"]
    centrality <- solve_block(I - G %*% (phi * I + gamma *
                                           e[["G_heterogeneity"]])) %*% e[["Ones"]]
  } else if (hypothesis == "het_l") {
    theta_0 <- coef(second_step)[names(coef(second_step)) %in% "theta_0"]
    theta_1 <- coef(second_step)[names(coef(second_step)) %in% "theta_1"]
    centrality <- solve_block(I - (theta_0 * I - theta_1 *
                                     e[["G_heterogeneity"]]) %*% G) %*% e[["Ones"]]
  } else if (hypothesis == "het_r") {
    eta_0 <- coef(second_step)[names(coef(second_step)) %in% "eta_0"]
    eta_1 <- coef(second_step)[names(coef(second_step)) %in% "eta_1"]
    centrality <- solve_block(I - G %*% (eta_0 * I - eta_1 *
                                           e[["G_heterogeneity"]])) %*% e[["Ones"]]

  } else if (hypothesis == "par") {
    phi_within <- coef(second_step)[names(coef(second_step)) %in% "phi_within"]
    phi_between <- coef(second_step)[names(coef(second_step)) %in% "phi_between"]
    centrality <- solve_block(I - phi_within * e[["G_within"]] -
                                phi_between * e[["G_between"]]) %*% e[["Ones"]]
  }

  return(centrality)

}
