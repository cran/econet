#' @keywords internal
#' @import minpack.lm
#' @noRd
try_boot_weight <- function(specification, sv, db, w) {

out <- tryCatch(
  {

    nlsLM(formula = specification, start = sv, data = db, trace= T, weights = w,
          control = nls.lm.control(ftol = sqrt(.Machine$double.eps) / 1000,
                                   ptol = sqrt(.Machine$double.eps)/1000,
                                   gtol = 0, diag = list(), epsfcn = 0,
                                   factor = 100, maxfev = integer(),
                                   maxiter = 300, nprint = 0))
  },
  error = function(cond){
    n <- length(sv)
    return(setNames(rep(NA,n),names(sv)))
  }
)
return(out)
}

