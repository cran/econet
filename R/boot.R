#' Bootstrap residuals with cross-sectional dependence
#' @param object numeric number
#' @param ... other arguments
#' @noRd
#' @export
"boot" <- function(object, ...){
  UseMethod("boot")
}
#' @rdname boot
#' @name boot
#' @title boot: Bootstrap residuals with cross-sectional dependence
#' @method boot econet
#' @aliases boot
#' @param object an object of class \code{econet} obtained using 'NLLS' estimation.
#' @param hypothesis string. One of \code{c("lim","het", "het_l", "het_r", "par")}.
#' @param group \code{NULL} or vector of positive integers specifying the indices for resampling within groups.
#' @param niter number of required iterations.
#' @param weights logical. It is \code{TRUE} if the object \code{object} was estimated using \code{weights}, and \code{FALSE} otherwise.
#' @param ... additional parameters
#' @return a numeric vector containing bootstrapped standard errors  (see Anselin, 1990). If the procedure is not feasible, it returns a vector of NAs.
#' @details
#' For additional details, see the vignette
#' Warning: This function is available only when \code{net_dep} is run with \code{estimation == "NLLS"}
#' @references Anselin, L., 1990, "Some robust approach to testing and estimation in spatial econometrics", Regional Science and Urban Economics, 20, 141-163.
#' @seealso \code{\link{net_dep}}
#' @examples
#' \donttest{
#' # Load data
#' data("db_cosponsor")
#' data("G_alumni_111")
#' db_model_B <- db_cosponsor
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
#' # object Linear-in-means model
#' lim_model_B <- net_dep(formula = f_model_B, data = db_model_B,
#'                        G = G_model_B, model = "model_B", estimation = "NLLS",
#'                        hypothesis = "lim", endogeneity = TRUE, correction = "heckman",
#'                        first_step = "standard",
#'                        exclusion_restriction = G_exclusion_restriction,
#'                        start.val = starting)
#' # Bootstrap
#' # Warning: this may take a very long time to run.
#' # Decrease the number of iterations to reduce runtime.
#' boot_lim_estimate <- boot(object = lim_model_B, hypothesis = "lim",
#'                           group = NULL, niter = 1000, weights = FALSE)
#' boot_lim_estimate
#' }
#' @importFrom "stats" "coef"  "constrOptim" "formula" "lm" "model.frame" "model.matrix" "model.response" "nlm" "nlminb" "optimize" "predict" "residuals" "sd" "setNames" "coefficients" "nobs" "pt"
#' @importFrom "utils" "data" "flush.console"
#' @export
"boot.econet" <- function(object, hypothesis = c("lim", "het", "het_l", "het_r", "par"),
                 group = NULL, niter, weights, ...) {

  if (class(object) == "econet") {
  object <- object[[1]]
  if (class(object) == "nls") {
    res<- list()
    for (i in 1:niter) {

      d <- sample_data(object, hypothesis, group)
      specification <- d[[1]]
      db <- d[[2]]
      sv <- replace_start_value(object, coef(object))

      if (weights == TRUE) {
        db$w <- object$data[["to_weight"]]
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
    res <- data.frame(coefficient = coefficients(object),
                      boot.Std.Error = res)
    res[, "boot.t.value"] <- res[, "coefficient"] / res[, "boot.Std.Error"]
    res[, "boot.p.value"] <- 2 * pt( - abs(res[, "boot.t.value"]),
                                     df = nobs(object) - nrow(res))

    return(res)

  } else {

    warning("The object was not estimated using NLLS")

  }
  }

  else warning(paste("This not an econet object or it was not estimated using 'NLLS'"))

}
