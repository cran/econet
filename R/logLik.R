#' logLik
#' @param object an object of class econet
#' @param centrality string. It is used when \code{object} is produced by \code{horse_race}.
#' @method logLik econet
#' @importFrom stats logLik
#' @noRd
#' @export
"logLik.econet" <- function(object, centrality = "parameter.dependent", ...) {

  if (!is.null(attributes(object)$attr)) {

    x <- object[[1]][[centrality]]

  } else {

    x <- object[[1]]

  }

  if (inherits(x, "nls") |
      inherits(x, "lm")) {
    res <- logLik(x)
  } else {
    res <- bbmle::logLik(x)
  }

  return(res)
}
#' @rdname logLik
#' @param object an object of class econet
#' @param centrality string. It is used when \code{object} is produced by \code{horse_race}.
#' @method logLik summary.econet
#' @importFrom stats logLik
#' @noRd
#' @export
"logLik.summary.econet" <- function(object, centrality = "parameter.dependent", ...) {

  if (!is.null(attributes(object)$attr)) {

    x <- object[[1]][[centrality]]

  } else {

    x <- object[[1]]

  }

  if (inherits(x, "nls") |
      inherits(x, "lm")) {
    res <- logLik(x)
  } else {
    res <- bbmle::logLik(x)
  }

  return(res)
}
