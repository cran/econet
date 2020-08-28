#' mle_party
#' @keywords internal
#' @noRd
#' @importFrom  bbmle parnames mle2
mle_party_split_with <- function(y, X, G_within_0, G_within_1, G_between,
                                 starting.values, boundL, boundU) {

  theta <- starting.values
  theta_L <- boundL
  theta_U <- boundU

  bbmle::parnames(ll_party_split_with) <- names(theta)

  fit <- mle2(ll_party_split_with, method="L-BFGS-B", start = theta,
              parnames = names(theta), vecpar = TRUE, lower = theta_L,
              upper = theta_U, control = list(maxit=1000),
              data = list(Y = y, X = X, G_between = G_between,
                          G_within_0 = G_within_0, G_within_1 = G_within_1))

  return(fit)
}

