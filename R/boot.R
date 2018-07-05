#' Bootstrap residuals with cross-sectional dependence
#'
#' @param fit First object in the list of outcomes returned by \code{net_dep}.
#' @param hypothesis string. One of \code{c("lim","het", "het_l", "het_r", "par")}.
#' @param group \code{NULL} or vector of positive integers specifying the indices for resampling within groups.
#' @param niter number of required iterations.
#' @param weights logical. It is \code{TRUE} if the object \code{fit} was estimated using \code{weights}, and \code{FALSE} otherwise.
#' @return a numeric vector containing bootstrapped standard errors  (see Anselin, 1990). If the procedure is not feasible, it returns a vector of NAs.
#' @details
#' For additional details, see the vignette
#' @references Anselin, L., 1990, "Some robust approach to testing and estimation in spatial econometrics", Regional Science and Urban Economics, 20, 141-163.
#' @seealso \code{\link{net_dep}}
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
#' # Bootstrap
#' # Warning: this may take a very long time to run.
#' # Decrease the number of iterations to reduce runtime.
#' boot_lim_estimate <- boot(fit = lim_estimate_model_B, hypothesis = "lim",
#'                           group = NULL, niter = 1000, weights = FALSE)
#' boot_lim_estimate
#' }
#' @importFrom "stats" "coef"  "constrOptim" "formula" "lm" "model.frame" "model.matrix" "model.response" "nlm" "nlminb" "optimize" "predict" "residuals" "sd" "setNames"
#' @importFrom "utils" "data" "flush.console"
#' @export
boot <- function(fit, hypothesis = c("lim", "het", "het_l", "het_r", "par"),
                 group = NULL, niter, weights) {

  res<- list()
  for (i in 1:niter) {

    d <- sample_data(fit, hypothesis, group)
    specification <- d[[1]]
    db <- d[[2]]
    sv <- replace_start_value(fit, coef(fit))

    if (weights == TRUE) {
      db$w <- fit$data[["to_weight"]]
    }

    if (weights == TRUE) {

      tmp <- try_boot_weight(specification, sv, db, w)

    } else {

      tmp <- try_boot(specification, sv, db)

    }

    if(all(is.na(tmp))){
      res[[i]] <- tmp
    } else {
      res[[i]] <- coef(tmp)
    }

  }

  res <- Reduce("cbind", res)
  res <- apply(res, 1, sd)
  return(res)
}
