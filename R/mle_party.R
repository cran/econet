#' mle_party
#' @keywords internal
#' @noRd
#' @importFrom  bbmle parnames mle2
mle_party <- function(y, X, G_within, G_between, starting.values, boundL,
                      boundU) {

  theta <- starting.values
  theta_L <- boundL
  theta_U <- boundU

  parnames(ll_party) <- names(theta)

  fit <- mle2(ll_party, method="L-BFGS-B", start = theta,
              parnames = names(theta), vecpar = TRUE, lower = theta_L,
              upper = theta_U, control = list(maxit=1000),
              data = list(Y = y, X = X, G_between = G_between, G_within = G_within))

  return(fit)
}

