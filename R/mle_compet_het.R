#' mle_compet_het
#' @keywords internal
#' @noRd
#' @importFrom  bbmle parnames mle2
mle_compet_het <- function(y, X, G, z, starting.values, boundL, boundU) {

  theta <- starting.values
  theta_L <- boundL
  theta_U <- boundU

  parnames(ll_compet_het_left) <- names(theta)
  X = X[,-1]

  fit <- mle2(ll_compet_het_left,  method="L-BFGS-B", start = theta,
              parnames = names(theta), vecpar = TRUE, lower = theta_L,
              upper = theta_U, control = list(maxit = 1000),
              data = list(Y = y, X = X, G = G, z = z))

  return(fit)

}
