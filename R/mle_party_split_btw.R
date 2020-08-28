#' mle_party
#' @keywords internal
#' @noRd
#' @importFrom  bbmle parnames mle2
mle_party_split_btw <- function(y, X, G_within,
                                G_between_01, G_between_10,
                                starting.values, boundL, boundU) {

  theta <- starting.values
  theta_L <- boundL
  theta_U <- boundU

  bbmle::parnames(ll_party_split_btw) <- names(theta)

  fit <- mle2(ll_party_split_btw, method="L-BFGS-B", start = theta,
              parnames = names(theta), vecpar = TRUE, lower = theta_L,
              upper = theta_U, control = list(maxit=1000),
              data = list(Y = y, X = X,
                          G_between_01 = G_between_01, G_between_10 = G_between_10,
                          G_within = G_within))

  return(fit)
}

