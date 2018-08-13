#' AIC
#' @param object an object of class econet
#' @param centrality string. It is used when \code{object} is produced by \code{horse_race}.
#' @method AIC econet
#' @importFrom stats AIC
#' @noRd
#' @export
"AIC.econet" <- function(object, centrality = "parameter.dependent", ..., k) {

  if (!is.null(attributes(object)$attr)) {
    x <- object[[1]][[centrality]]
  } else {
    x <- object[[1]]
  }

  res <- AIC(x)

  return(res)
}
#' @rdname AIC
#' @param object an object of class summary.econet
#' @param centrality string. It is used when \code{object} is produced by \code{horse_race}.
#' @method AIC summary.econet
#' @importFrom stats AIC
#' @noRd
#' @export
"AIC.summary.econet" <- function(object, centrality = "parameter.dependent", ..., k) {

  if (!is.null(attributes(object)$attr)) {
    x <- object[[1]][[centrality]]
  } else {
    x <- object[[1]]
  }

  res <- AIC(x)

  return(res)
}
