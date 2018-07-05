#' Quantification of marginal effects in linear-in-means models.
#'
#' @param fit first object in the list of outcomes returned by \code{net_dep} (available only if the argument \code{model} is set to \code{"model_B"}).
#' @param Gn an object of class \code{Matrix} representing the social network.
#' @return an object of class \code{data.frame} listing direct and indirect variable effects (mean, standard deviation, max, min).
#' @details
#' For additional details, see the vignette
#' @examples
#' \donttest{
#' # Load data
#' data("db_cosponsor")
#' data("G_alumni_111")
#' db_model_B <- subset(db_cosponsor, time == 3)
#' G_model_B <- G_cosponsor_111
#' G_exclusion_restriction <- G_alumni_111
#' are_factors <- c("party", "gender", "nchair")
#' db_model_B[are_factors] <- lapply(db_model_B[are_factors], factor)
#'
#' # Specify formula
#' f_model_B <- formula("les ~gender + party + nchair")
#'
#' # Specify starting values
#' starting <- c(alpha = 0.23952,
#'               beta_gender1 = -0.22024,
#'               beta_party1 = 0.42947,
#'               beta_nchair1 = 3.09615,
#'               phi = 0.40038,
#'               unobservables = 0.07714)
#'
#' # Fit Linear-in-means model
#' lim_model_B <- net_dep(formula = f_model_B, data = db_model_B,
#'                        G = G_model_B, model = "model_B", estimation = "NLLS",
#'                        hypothesis = "lim", endogeneity = TRUE, correction = "heckman",
#'                        first_step = "standard",
#'                        exclusion_restriction = G_exclusion_restriction,
#'                        start.val = starting)
#' # Store and print results
#' lim_estimate_model_B <- lim_model_B[[1]]; summary(lim_estimate_model_B)
#' # Quantification
#' quantify(fit = lim_estimate_model_B, Gn = G_model_B)
#' }
#' # WARNING, This toy example is provided only for runtime execution.
#' # Please refer to previous examples for sensible calculations.
#' data("db_alumni_test")
#' data("G_model_A_test")
#' db_model_A <- db_alumni_test
#' G_model_A <- G_model_A_test
#' f_model_A <- formula("les ~ dw")
#' lim_model_A_test <- net_dep(formula = f_model_A, data = db_model_A,
#'                        G = G_model_A, model = "model_A", estimation = "NLLS",
#'                        hypothesis = "lim", start.val = c(alpha = 0.09030594,
#'                                                          beta_dw = 1.21401940,
#'                                                          phi = 1.47140647))
#' quantify(lim_model_A_test[[1]], G_model_A)
#' @seealso \code{\link{net_dep}}
#' @importFrom Matrix solve
#' @export
quantify <- function(fit, Gn) {

  if (class(fit) == "mle2") {
    fit <- summary(fit)@coef
  } else if (class(fit) == "nls") {
    fit <- summary(fit)$coefficients
  }

  fit <- fit[ - which(rownames(fit)%in%"alpha"), ]
  sel <- which(rownames(fit) %in% "phi")
  regressor <- fit[ - sel, 1]
  parameter <- fit[sel, 1]
  to_display <- rownames(regressor)
  fit <- c(parameter,regressor)

  res.b <- fit
  n <- nrow(Gn)

  D <- NULL
  M <- NULL
  Sn <- solve(diag(dim(Gn)[1]) - res.b[1] * Gn) %*% diag(n)

  for (i in 2:length(res.b)) {
    beta_k <- i

    Me <- rep(0, n)
    Ie <- rep(0, n)
    Dx <- rep(1, n)
    C_i <- Sn * res.b[beta_k]
    Me <- diag(C_i)
    D <- rbind(D, matrix(c(mean(Me), sd(Me), max(Me), min(Me)), 1, 4))

    Me2 <- Me
    Me <- diag(length(Me2))
    diag(Me) <- Me2

    C_i <- C_i - Me;
    Ie <- c(C_i[C_i != 0])
    M <- rbind(M, matrix(c(mean(Ie), sd(Ie), max(Ie), min(Ie)), 1, 4))
  }

  res <- cbind(round(D,4), round(M,4))
  res <- data.frame(res)
  res <- cbind(round(regressor,4), res)
  colnames(res) <- c("beta",
                     paste0(rep("Direct_"),
                            c("mean", "std", "max", "min")),
                     paste0(rep("Indirect_"),
                            c("mean", "std", "max", "min")))
  rownames(res) <- to_display

  return(res)
}




