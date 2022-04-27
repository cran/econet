#' formula
#' @param x an object of class econet
#' @param print string. If \code{"main.equation"} returns estimates for the main equation. \cr
#' If \code{"first.step"} returns estimates for the first step.
#' If \code{"centrality"} returns quantile distribution for the estimated centrality.
#' @param centrality string. It is used when \code{object} is produced by \code{horse_race}.
#' @method formula econet
#' @noRd
#' @export
"formula.econet" <- function(x, print = "main.equation", centrality = "parameter.dependent", ...) {

  if (print == "main.equation") {

    if (!is.null(attributes(x)$attr)) {

      x <- x[[1]][[centrality]]

    } else {

      x <- x[[1]]

    }

    if (inherits(x, "nls")) {
      if (!is.null(x$m$mformula)) {
        res <- formula(x$m$mformula)
      } else {
        res <- formula(x, ...)
      }
    } else {
      res <- bbmle::formula(x, ...)
    }

  } else if(print == "first.step") {
    res <- formula(x[[2]], ...)
  }

  class(res) <- "formula.econet"
  return(res)
}
