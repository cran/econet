#' Maximum likelihood function for model het_l
#' @keywords internal
#' @noRd
#' @importFrom utils flush.console
ll_het_left <- function(params){

  name.params <- names(params)
  beta <- params[ - which(name.params %in% c("theta_0","theta_1","sigma"))]
  beta <- unlist(beta)
  sigma <- unlist(params[which(name.params %in% "sigma")])
  theta_0 <- unlist(params[which(name.params %in% "theta_0")])
  theta_1 <- unlist(params[which(name.params %in% "theta_1")])

  cat("pars: ", theta_0, theta_1, sigma, beta, '\n')
  flush.console()

  n <- nrow(G)
  I <- diag(n)
  X <- as.matrix(X)
  z <- diag(z)

  sigma_2 <- sigma
  reg_par <- solve_block(I - G %*% (theta_0 * I + theta_1 * z))
  reg_par_2 <- solve_block(I - (theta_0 * I + theta_1 * z) %*% t(G))
  omega <- sigma_2 * reg_par %*% reg_par_2

  det_reg_par <- as.numeric(determinant(reg_par %*% reg_par_2,
                                        logarithm = T)$modulus)

  other_params <- Y - reg_par %*% (X %*% beta)

  res <- -(n / 2) * log(sigma_2) - (1 / 2) * det_reg_par - (1 / 2) *
    t(other_params) %*% solve_block(omega) %*% other_params

  print( - res)
}
