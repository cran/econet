#' Maximum likelihood function for model het_r
#' @keywords internal
#' @noRd
#' @importFrom utils flush.console
ll_het_right <- function(params){

  name.params <- names(params)
  beta <- params[ - which(name.params %in% c("eta_0", "eta_1", "sigma"))]
  beta <- unlist(beta)
  sigma <- unlist(params[which(name.params %in% "sigma")])
  eta_0 <- unlist(params[which(name.params %in% "eta_0")])
  eta_1 <- unlist(params[which(name.params %in% "eta_1")])

  cat("pars: ", eta_0, eta_1, sigma,  beta, '\n')
  flush.console()

  n <- nrow(G)
  I <- diag(n)
  X <- as.matrix(X)
  z <- diag(z)

  sigma_2 <- sigma
  reg_par <- solve_block(I - (eta_0 * I + eta_1 * z) %*% G)
  reg_par_2 <- solve_block(I - (eta_0 * I + eta_1 * z) %*% t(G))
  omega <- sigma_2 * reg_par %*% reg_par_2

  det_reg_par <- as.numeric(determinant(reg_par %*% reg_par_2,
                                        logarithm = T)$modulus)

  other_params <- Y - reg_par %*% (X %*% beta)

  res <- - (n / 2) * log(sigma_2) - (1 / 2) * det_reg_par - (1 / 2) *
    t(other_params) %*% solve_block(omega) %*% other_params

  print( - res)
}
