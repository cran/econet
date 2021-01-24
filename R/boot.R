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
#' @param hypothesis string. One of \code{c("lim","het", "het_l", "het_r", "par", "par_split_with", "par_split_btw", "par_split_with_btw")}.
#' @param group \code{NULL} or vector of positive integers specifying the indices for resampling within groups.
#' @param niter number of required iterations.
#' @param weights logical. It is \code{TRUE} if the object \code{object} was estimated using \code{weights}, and \code{FALSE} otherwise. Default \code{FALSE}.
#' @param na.rm logical. Should missing values (including \code{NaN}) be removed?
#' @param delta Default is \code{NULL}.It has to be a number between zero (included) and one (excluded). When used, \code{econet} performs a constrained NLLS estimation. In this case, the estimated peer effect parameter, taken in absolute value, is forced to be higher than zero and lower than the spectral radius of \code{G}. Specifically, \code{delta} is a penalizing factor, decreasing the goodness of fit of the NLLS estimation, when the peer effect parameter approaches one of the two bounds. Observe that very high values of \code{delta} may cause NLLS estimation not to converge.
#' @param parallel logical. It is \code{TRUE} if the user wants to make use of parallelization. Default is \code{FALSE}. Observe that this option is still in its beta version and it has been tested only for the Windows platform.
#' @param cl numeric. Number of cores to be used for parallelization.
#' @param ... additional parameters
#' @return a numeric vector containing bootstrapped standard errors  (see Anselin, 1990). If the procedure is not feasible, it returns a vector of NAs.
#' @details
#' For additional details, see the vignette.
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
#' # If you run econet on a Windows platform, you can try to set the
#' # argument parallel = TRUE. However note that this option is still
#' # in its beta version.
#' boot_lim_estimate <- boot(object = lim_model_B, hypothesis = "lim",
#'                           group = NULL, niter = 10, weights = FALSE)
#' boot_lim_estimate
#' }
#' @importFrom "stats" "coef"  "constrOptim" "formula" "lm" "model.frame" "model.matrix" "model.response" "nlm" "nlminb" "optimize" "predict" "residuals" "sd" "setNames" "coefficients" "nobs" "pt"
#' @importFrom "utils" "data" "flush.console" "txtProgressBar" "setTxtProgressBar"
#' @importFrom "doParallel" "registerDoParallel"
#' @import "parallel"
#' @import "foreach"
#' @import "progressr"
#' @export
"boot.econet" <- function(object, hypothesis = c("lim", "het", "het_l", "het_r", "par", "par_split_with", "par_split_btw", "par_split_with_btw"),
                 group = NULL, niter, weights = FALSE, delta = NULL, na.rm = FALSE, parallel = FALSE, cl, ...) {

  if (class(object) == "econet") {
  object <- object[[1]]
  if (class(object) == "nls") {

    if (parallel == FALSE & is.null(delta)) {
    res<- list()
    for (i in 1:niter) {

      d <- sample_data(object, hypothesis, group, delta = object$delta)
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

    } else if (parallel == FALSE & !is.null(delta)) {
      res<- list()
      for (i in 1:niter) {

        d <- sample_data(object, hypothesis, group, delta = object$delta)
        specification <- d[[1]]
        db <- d[[2]]
        sv <- replace_start_value(object, coef(object))

        if (weights == TRUE) {
          nls_const <- nls_constrained(nls_formula_fit = specification, G = db$G, e = db,
                                       hypothesis = hypothesis, model = "model_B", delta = delta, to_weight = db$w)

          specification <- nls_const$nls_formula_fit
          db <- nls_const$e
          db$w <- db[["to_weight"]]

          tmp <- try_boot_weight(specification, sv, db, w)

        } else {

          nls_const <- nls_constrained(nls_formula_fit = specification, G = db$G, e = db,
                                       hypothesis = hypothesis, model = "model_B", delta = delta, to_weight = db$w)

          specification <- nls_const$nls_formula_fit
          db <- nls_const$e
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

    } else if (parallel == TRUE & is.null(delta)) {
      cl <- makeCluster(cl)
      doParallel::registerDoParallel(cl)
      iterations <- niter
      progressr::with_progress({
        p <- progressr::progressor(along = iterations)

        res <- foreach(i = 1:iterations, .export = c("select_X", "select_beta",
                                                     "par_dep", "dep_mat", "par_dep_mat",
                                                     "sample_residuals", "create_y_hat",
                                                     "sample_data", "replace_start_value",
                                                     "try_boot_weight", "try_boot", "solve_block"),
                       .packages = c("igraph", "MASS", "minpack.lm")) %dopar%
          {
            d <- sample_data(object, hypothesis, group, delta = object$delta)
            specification <- d[[1]]
            db <- d[[2]]
            sv <- replace_start_value(object, coef(object))

            if (weights == TRUE) {
              db$w <- object$data[["to_weight"]]
            }

            if (weights == TRUE) {

              tmp <- try_boot_weight(specification, sv, db, db$w)

            } else {

              tmp <- try_boot(specification, sv, db)

            }

            if(all(is.na(tmp))){
              res <- tmp
            } else {
              res <- coef(tmp)
            }

            return(res)
          }
      })


      res <- Reduce("cbind", res)
      res <- apply(res, 1, sd, na.rm = na.rm)
      res <- data.frame(coefficient = coefficients(object),
                        boot.Std.Error = res)
      res[, "boot.t.value"] <- res[, "coefficient"] / res[, "boot.Std.Error"]
      res[, "boot.p.value"] <- 2 * pt( - abs(res[, "boot.t.value"]),
                                       df = nobs(object) - nrow(res))

      stopCluster(cl)
    } else if (parallel == TRUE & !is.null(delta)) {
      cl <- makeCluster(cl)
      doParallel::registerDoParallel(cl)
      iterations <- niter
      progressr::with_progress({
        p <- progressr::progressor(along = iterations)

        res <- foreach(i = 1:iterations, .export = c("select_X", "select_beta",
                                                     "par_dep", "dep_mat", "par_dep_mat",
                                                     "sample_residuals", "create_y_hat",
                                                     "sample_data", "replace_start_value",
                                                     "try_boot_weight", "try_boot", "solve_block",
                                                     "nls_constrained"),
                       .packages = c("igraph", "MASS", "minpack.lm")) %dopar%
          {
            d <- sample_data(object, hypothesis, group, delta = object$delta)
            specification <- d[[1]]
            db <- d[[2]]
            sv <- replace_start_value(object, coef(object))

            if (weights == TRUE) {

              nls_const <- nls_constrained(nls_formula_fit = specification, G = db$G, e = db,
                                           hypothesis = hypothesis, model = "model_B", delta = delta, to_weight = db$w)

              specification <- nls_const$nls_formula_fit
              db <- nls_const$e
              db$w <- db[["to_weight"]]

              tmp <- try_boot_weight(specification, sv, db, db$w)

            } else {

              nls_const <- nls_constrained(nls_formula_fit = specification, G = db$G, e = db,
                                           hypothesis = hypothesis, model = "model_B", delta = delta, to_weight = db$w)

              specification <- nls_const$nls_formula_fit
              db <- nls_const$e

              tmp <- try_boot(specification, sv, db)

            }

            if(all(is.na(tmp))){
              res <- tmp
            } else {
              res <- coef(tmp)
            }

            return(res)
          }
      })


      res <- Reduce("cbind", res)
      res <- apply(res, 1, sd, na.rm = na.rm)
      res <- data.frame(coefficient = coefficients(object),
                        boot.Std.Error = res)
      res[, "boot.t.value"] <- res[, "coefficient"] / res[, "boot.Std.Error"]
      res[, "boot.p.value"] <- 2 * pt( - abs(res[, "boot.t.value"]),
                                       df = nobs(object) - nrow(res))

      stopCluster(cl)
    }

    return(res)

  } else {

    warning("The object was not estimated using NLLS")

  }
  }

  else warning(paste("This not an econet object or it was not estimated using 'NLLS'"))

}
