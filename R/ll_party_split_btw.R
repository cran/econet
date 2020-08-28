#' Maximum likelihood function for model party
#' @keywords internal
#' @noRd
#' @importFrom utils flush.console
ll_party_split_btw <- function(params){

  name.params <- names(params)
  beta <- params[ - which(name.params %in% c("phi_within", "phi_between_01", "phi_between_10",
                                             "sigma"))]
  beta <- unlist(beta)
  sigma <- unlist(params[which(name.params %in% "sigma")])
  phi_within <- unlist(params[which(name.params %in% "phi_within")])
  phi_between_01 <- unlist(params[which(name.params %in% "phi_between_01")])
  phi_between_10 <- unlist(params[which(name.params %in% "phi_between_10")])

  cat("pars: ", phi_within, phi_between_01, phi_between_10, sigma, beta, '\n')
  flush.console()

  n <- nrow(G_within)
  I <- diag(n)
  X <- as.matrix(X)

  sigma_2 <- sigma
  reg_par <- solve_block(I - (phi_within * G_within) - (phi_between_01 * G_between_01) - (phi_between_10 * G_between_10))
  reg_par_2 <- solve_block(I - (phi_within * t(G_within)) -
                             (phi_between_01 * t(G_between_01)) -
                             (phi_between_10 * t(G_between_10)))
  omega <- sigma_2 * reg_par %*% reg_par_2
  omega_inv = (I - (phi_within * t(G_within)) - (phi_between_01 * t(G_between_01)) - (phi_between_10 * t(G_between_10))) *
    (I - (phi_within * G_within) - (phi_between_01 * G_between_01) - (phi_between_10 * G_between_10)) / sigma_2

  det_reg_par <- as.numeric(determinant(reg_par %*% reg_par_2,
                                        logarithm = T)$modulus)

  other_params <- Y - reg_par %*% (X %*% beta)

  res <- - (n / 2) * log(sigma_2) - (1 / 2) * det_reg_par - (1/2) *
    t(other_params) %*% omega_inv %*% other_params

  print( - res)
}
