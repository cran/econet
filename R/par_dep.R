#' Extract estimated decay-parameter
#'
#' @param fit First object in the list of outcomes returned by net_dep.
#' @param type hypothesis: e.g. \code{"lim"}, \code{"het"}, \code{"het_l"}, \code{"het_r"}, \code{"par"}.
#' @return Value of the estimated decay-parameter.
#' @keywords internal
#' @importFrom stats coef
#' @noRd
par_dep <- function(fit, type){
  switch(type,
         lim = coef(fit)[names(coef(fit)) %in% "phi"],
         het = coef(fit)[names(coef(fit)) %in% c("phi", "gamma")],
         het_l = coef(fit)[names(coef(fit)) %in% c("theta_0", "theta_1")],
         het_r = coef(fit)[names(coef(fit)) %in% c("eta_0", "eta_1")],
         par = coef(fit)[names(coef(fit)) %in% c("phi_within", "phi_between")],
         par_split_btw = coef(fit)[names(coef(fit)) %in% c("phi_within", "phi_between_01", "phi_between_10")],
         par_split_with = coef(fit)[names(coef(fit)) %in% c("phi_within_0", "phi_within_1", "phi_between")],
         par_split_with_btw = coef(fit)[names(coef(fit)) %in% c("phi_within_0", "phi_within_1", "phi_between_01", "phi_between_10")]
  )
}
