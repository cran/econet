#' print
#' @aliases print.econet
#' @method print econet
#' @param x an object of class econet
#' @param digits minimum number of significant digits to be used for most numbers.
#' @noRd
#' @rdname print
#' @export
"print.econet" <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {

  main_equation <- x[[1]]

  if( inherits(main_equation, "nls") |
      inherits(main_equation, "lm") ) {

    print.default(format(coef(main_equation), digits = digits), print.gap = 2L,
                  quote = FALSE)
    cat("\n")
    invisible(main_equation)
  } else {

    print.default(format(bbmle::coef(main_equation), digits = digits), print.gap = 2L,
                  quote = FALSE)
    cat("\n")
    invisible(main_equation)
  }
}
#' @rdname print
#' @param x numeric number
#' @param ... other arguments
#' @method print summary.econet
#' @noRd
#' @export
"print.summary" <- function(x, ...){
  UseMethod("print.summary")
}
#' @rdname print
#' @param x an object of class econet
#' @param print string. If \code{"main.equation"} returns estimates for the main equation. \cr
#' If \code{"first.step"} returns estimates for the first step.
#' If \code{"centrality"} returns quantile distribution for the estimated centrality.
#' @param centrality string. It is used when \code{object} is produced by \code{horse_race}.
#' @param digits minimum number of significant digits to be used for most numbers.
#' @param signif.stars logical; if \code{TRUE}, P-values are additionally encoded visually as ‘significance stars’ in order to help scanning of long coefficient tables. It defaults to the \code{show.signif.stars} slot of \code{options}.
#' @noRd
#' @method print summary.econet
#' @export
"print.summary.econet" <- function(x, print = "main.equation",
                                   centrality = "parameter-dependent",
                                   digits = max(3L, getOption("digits") - 3L),
                                   signif.stars = getOption("show.signif.stars"),
                                   ...) {


  if(is.null(x)) {
    warning("This  object was estimated with endogeneity = 'FALSE'")
  } else {
    if (attributes(x)$attr == "main.equation") {

      formula <- x$formula
      aic <- x$aic
      loglik <- x$loglik
      coefficients <- x$coefficients

    } else if (attributes(x)$attr == "first.step") {

      formula <- x$formula
      coefficients <- x$coefficients
      r2 <- x$r2

    } else if (attributes(x)$attr == "centrality") {

      res <- x$centrality.dist

    }

    if (attributes(x)$attr == "main.equation") {
      cat("Call:\n")
      cat(paste("Main Equation: ", paste(formula[2], formula[3], sep = ' ~ ')))
      cat("\n")

      if( inherits(x, "nls") |
          inherits(x, "lm") ) {

        printCoefmat(coefficients, digits = digits, signif.stars = signif.stars,
                     na.print = "NA", ...)

      } else {

        printCoefmat(coefficients, digits = digits, signif.stars = signif.stars,
                     na.print = "NA", ...)

      }

      cat("\n")
      cat(paste0("AIC: ", round(aic, 2)),
          paste0(" loglik: ", round(loglik, 2)))
      invisible(x)

    } else if (attributes(x)$attr == "first.step" & !is.null(x)) {

      cat(paste("First step: ", paste(formula[2], formula[3], sep = ' ~ ')))
      cat("\n")
      printCoefmat(coefficients, digits = digits, signif.stars = signif.stars,
                   na.print = "NA", ...)
      cat("\n")
      cat(paste0("R2: ", round(r2, 2)))
      cat("\n")
      invisible(x)

    } else if (attributes(x)$attr == "centrality") {

      cat("Call:\n")
      cat(paste0(x$centrality, " centrality distribution"))
      cat("\n")
      print.default(format(res), print.gap = 2L,
                    quote = FALSE)
      cat("\n")
      invisible(x)

    }

  }

}
