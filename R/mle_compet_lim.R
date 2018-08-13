#' mle_compet_lim
#' @keywords internal
#' @noRd
#' @importFrom  bbmle parnames mle2
mle_compet_lim <- function(y, X, G, z, starting.values, boundL, boundU) {

  theta <- starting.values
  theta_L <- boundL
  theta_U <- boundU
  bbmle::parnames(ll_compet_lim) <- names(theta)
  X = X[, -1]
  fit <- mle2(ll_compet_lim, method = "L-BFGS-B", start = theta,
              parnames = names(theta), vecpar = TRUE, lower = theta_L,
              upper = theta_U, control = list(maxit = 1000),
              data = list(Y = y, X = X, G = G, z = z))

  return(fit)

}
