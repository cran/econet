#' vcov
#' @param print string. If \code{"main.equation"} returns estimates for the main equation. \cr
#' If \code{"first.step"} returns estimates for the first step.
#' If \code{"centrality"} returns quantile distribution for the estimated centrality.
#' @param centrality string. It is used when \code{object} is produced by \code{horse_race}.
#' @method vcov econet
#' @importFrom stats vcov
#' @noRd
#' @export
"vcov.econet" <- function(object, print = "main.equation", centrality = "parameter.dependent", ...) {

  x <- object

  if(print == "main.equation") {

    if (!is.null(attributes(x)$attr)) {

      x <- x[[1]][[centrality]]

    } else {

      x <- x[[1]]

    }

    if (inherits(x, "nls") |
        inherits(x, "lm")) {
      res <- vcov(x, ...)
    } else {
      res <- bbmle::vcov(x, ...)
    }

  } else if(print == "first.step") {
    res <- vcov(x[[2]], ...)
  }

  class(res) <- "vcov.econet"
  return(res)
}
