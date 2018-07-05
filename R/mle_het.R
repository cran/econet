#' mle_compet_lim
#' @keywords internal
#' @noRd
#' @importFrom  bbmle parnames mle2
mle_het <- function(y, X, G, z, side, starting.values, boundL, boundU) {

  theta <- starting.values
  theta_L <- boundL
  theta_U <- boundU

  if(side == "het_l"){
    parnames(ll_het_left) <- names(theta)
    fit <- mle3(ll_het_left, optimizer = "constrOptim",
                method = "Nelder-Mead", start = theta, parnames = names(theta),
                ui = rbind(c(rep(0, length(theta) - 3), - 1, - 1,0),
                           c(rep(0, length(theta) - 3), 1, 1, 0),
                           c(rep(0, length(theta) - 1), 1)),
                ci = c( - 1 / max(colSums(G)), - 1 / max(colSums(G)), 0),
                vecpar = TRUE, skip.hessian = TRUE,
                data = list(Y = y, X = X, G = G, z = z))

   } else if (side == "het_r"){
    parnames(ll_het_right) <- names(theta)
    fit <- mle3(ll_het_right, optimizer = "constrOptim",
                method = "Nelder-Mead", start = theta, parnames = names(theta),
                ui = rbind(c(rep(0, length(theta) - 3), - 1, - 1, 0),
                           c(rep(0, length(theta) - 3), 1, 1, 0),
                           c(rep(0, length(theta) - 1), 1)),
                ci = c( - 1 / max(rowSums(G)), - 1 / max(rowSums(G)), 0),
                vecpar = TRUE, skip.hessian = TRUE,
                data = list(Y = y, X = X, G = G, z = z))
   }

  return(fit)
}
