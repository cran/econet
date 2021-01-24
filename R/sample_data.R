#' Generate new dependent variable by resampling model's predicted values
#'
#' @param fit First object in the list of outcomes returned by net_dep.
#' @param type hypothesis: e.g. \code{"lim"}, \code{"het"}, \code{"het_l"}, \code{"het_r"}, \code{"par"}.
#' @param group \code{NULL} or vector of positive integers specifying the indices for resampling within groups.
#' @param delta logical. It is set to \code{TRUE} if the object was estimated using a barrier function, and \code{FALSE} otherwise.
#' @return List of two objects: i) model specification, ii) new dataset.
#' @keywords internal
#' @importFrom stats coef
#' @noRd
sample_data <- function(fit, hypothesis, group, delta){

  if (delta) {
    specification <- fit$m$mformula
  } else {
    specification <- formula(fit)
  }

  y_name <- specification[[2]]
  y <- fit$data[[as.character(y_name)]]

  X_names <- select_X(specification)
  X <- Reduce("cbind",fit$data[names(fit$data) %in% X_names])
  if(is.null(dim(X))){ X <- matrix(X) }
  colnames(X) <- X_names

  beta_names <- select_beta(specification)
  beta <- coef(fit)
  beta <- beta[names(beta) %in% beta_names]
  beta <- beta[order(match(names(beta), beta_names))]

  res <- residuals(fit)

  if (delta) {
    res <- res[-length(res)]
  }

  hypothesis <- hypothesis
  pd <- par_dep(fit, hypothesis)
  mat <- dep_mat(fit, hypothesis)
  n <- length(res)
  I <- diag(n)
  pdm <- par_dep_mat(par_dep = pd, G = mat, I = I, type = hypothesis)

  sample_res <- sample_residuals(res = res, par_dep_mat = pdm, group = group)
  y_hat <- create_y_hat(beta = beta, X = X, par_dep_mat = pdm,
                        sample_res = sample_res)

  fit$data[[y_name]] <- y_hat

  list(formula = specification, data = fit$data)
}
