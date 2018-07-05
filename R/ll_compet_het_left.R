#' Maximum likelihood function for model compet_het
#' @keywords internal
#' @noRd
#' @importFrom utils flush.console
ll_compet_het_left <- function(params){

  name.params <- names(params)
  beta <- params[ - which(name.params %in% c("alpha", "phi", "gamma", "sigma"))]
  beta <- unlist(beta)
  sigma <- unlist(params[which(name.params %in% "sigma")])
  phi <- unlist(params[which(name.params%in%"phi")])
  gamma <- unlist(params[which(name.params%in%"gamma")])
  alpha <- unlist(params[which(name.params%in%"alpha")])

  cat("pars: ", alpha, phi, gamma, sigma, beta, '\n')
  flush.console()

  n <- nrow(G)
  I <- diag(n)
  X <- as.matrix(X)
  z <- diag(z)

  sigma_2 <- sigma
  reg_par <- solve_block(I - G %*% (phi * I + gamma * z))

  other_params <- Y - alpha * (reg_par %*% rep(1, n)) - (X %*% beta)

  res <- - (n / 2) * log(sigma_2) - (1 / 2) *
    t(other_params) %*% other_params / sigma_2

  print( - res)
}
