#' summary
#' @param object an object of class \code{econet}.
#' @param print string. If \code{"main.equation"} returns estimates for the main equation. \cr
#' If \code{"first.step"} returns estimates for the first step.
#' If \code{"centrality"} returns quantile distribution for the estimated centrality.
#' @param centrality string. It is used when \code{object} is produced by \code{horse_race}.
#' @param digits minimum number of significant digits to be used for most numbers.
#' @param signif.stars logical; if \code{TRUE}, P-values are additionally encoded visually as ‘significance stars’ in order to help scanning of long coefficient tables. It defaults to the \code{show.signif.stars} slot of \code{options}.
#' @method summary econet
#' @importFrom stats AIC printCoefmat vcov
#' @noRd
#' @export
"summary.econet" <- function(object, print = "main.equation",
                             centrality = "parameter.dependent",
                             digits = max(3L, getOption("digits") - 3L),
                             signif.stars = getOption("show.signif.stars"), ...) {

    if (print == "main.equation") {

      main_equation <- object[[1]]

      if (!is.null(attributes(object)$attr)) {
        main_equation <- main_equation[[centrality]]
      }

      formula <- formula(object, print, centrality, ...)
      aic <- AIC(object, centrality)
      loglik <- logLik(object, centrality)

      if( class(main_equation) == "nls" |
          class(main_equation) == "lm" ) {

        coefficients <- summary(main_equation)$coefficients

      } else {

        coefficients <- bbmle::summary(main_equation)@coef

      }

      res <- list(formula = formula,
                  aic = aic, loglik = loglik,
                  coefficients = coefficients,
                  print = print,
                  centrality = centrality,
                  centralities = object[[2]])

      class(res) <- "summary.econet"
      attributes(res)$attr <- "main.equation"

    } else if (print == "first.step") {

      if(!is.null(object[[3]])) {

        first_step <- object[[3]]

        if (!is.null(attributes(object)$attr)) {
          first_step <- first_step[[centrality]]
        }

        formula <- formula(first_step)
        r2 <- summary(first_step)$r.squared
        coefficients <- summary(first_step)$coefficients
        res <- list(formula = formula,
                    r2 = r2,
                    coefficients = coefficients,
                    print = print,
                    centrality = centrality,
                    centralities = object[[2]])

        class(res) <- "summary.econet"
        attributes(res)$attr <- "first.step"

      } else {

        res <- NULL

      }

    } else if (print == "centrality") {

      object <- object[[2]][, centrality]
      res <- list(centrality.dist = summary(c(object)), print = print, centrality = centrality)

      class(res) <- "summary.econet"
      attributes(res)$attr <- "centrality"

    }

  res

}


