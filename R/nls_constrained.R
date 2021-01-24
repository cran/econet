#' Create a list of objects to be passed to the function \code{nlsLM}
#' @param nls_formula_fit formula an object of class \code{formula}: a symbolic description of the model to be fitted. The constant (i.e. intercept) and the autogressive parameter needs not to be specified.
#' @param G G an object of class \code{Matrix} representing the social network.
#' @param e list of objects created either by the user (e.g. \code{G}) or by the function \code{net_dep} (e.g. \emph{"unobservables"})
#' @param model string. One of \code{c("model_A","model_B")}.
#' @param hypothesis string. It indicates the how peer effects should be modeled.
#' @param delta Default is \code{NULL}. To be used when \code{estimation = "NLLS"}. It has to be a number between zero (included) and one (excluded). When used, \code{econet} performs a constrained NLLS estimation. In this case, the estimated peer effect parameter, taken in absolute value, is forced to be higher than zero and lower than the spectral radius of \code{G}. Specifically, \code{delta} is a penalizing factor, decreasing the goodness of fit of the NLLS estimation, when the peer effect parameter approaches one of the two bounds. Observe that very high values of \code{delta} may cause NLLS estimation not to converge.
#' @param to_weight an optional vector of weights to be used in the fitting process to indicate that different observations have different variances. Should be \code{NULL} or a numeric vector. If non-\code{NULL}, it can be used to fit a weighted non-linear least squares (\code{estimation = "NLLS"}).
#' @return a list of two objects: i) an object of class \code{formula}; ii) a list of objects containing the data to be used in the regression.
#' @importFrom formula.tools rhs.vars lhs.vars
#' @keywords internal
#' @noRd
nls_constrained <- function(nls_formula_fit, G, e, hypothesis, model, delta, to_weight) {

  egv <- eigen(G)
  dim_barrier <- function(hypothesis) {
    switch(hypothesis,
           lim = 1,
           het = 2,
           het_l = 2,
           het_r = 2,
           par = 2,
           par_split_with = 3,
           par_split_btw = 3,
           par_split_with_btw = 4)
  }

  e_new <- e
  max_egv <- max(abs(egv$values))
  e_new$max_egv <- max_egv
  n <- nrow(e_new$I)
  e_new$n <- n
  e_new$I <- rbind(e_new$I, matrix(0, nrow = dim_barrier(hypothesis), ncol = n))
  e_new$I <- cbind(e_new$I, matrix(0, ncol = dim_barrier(hypothesis),
                                   nrow = nrow(e_new$I)))
  e_new$G <- rbind(e_new$G, matrix(0, nrow = dim_barrier(hypothesis),
                                   ncol = n))
  e_new$G <- cbind(e_new$G, matrix(0, ncol = dim_barrier(hypothesis),
                                   nrow = nrow(e_new$I)))
  e_new$delta <- delta

  sel_slot <- which(names(e_new) %in% c("max_egv", "delta", "n", "I", "G",
                                        "G_between", "G_within", "G_between_01",
                                        "G_between_10", "G_within_0", "G_within_1",
                                        "G_heterogeneity"))
  e_new[-sel_slot] <- lapply(e_new[-sel_slot], function(x) {
    if (is.numeric(x)) {
      x <- c(x, rep(0, dim_barrier(hypothesis)))
    } else if (is.factor(x)) {
      x <- c(x, rep("NA", dim_barrier(hypothesis)))
    }
    return(x)
  })

  if (!is.null(to_weight)) {
    new_to_weight <- c(to_weight, rep(to_weight[1], dim_barrier(hypothesis)))
  } else {
    new_to_weight <- NULL
  }

  if (model == "model_B") {
    rhs <- rhs.vars(nls_formula_fit)
    lin_rhs <- strsplit(rhs, "alpha")[[1]][2]
    lin_rhs <- paste0("(alpha", lin_rhs)
    lin_rhs <- gsub(" \\+ ", "[1:n] + ", lin_rhs)
    lin_rhs <- gsub("\\)", "[1:n])", lin_rhs)
    sp_rhs <- strsplit(rhs, "alpha")[[1]][1]
    sp_rhs <- gsub(" \\(", "", sp_rhs)
    sp_rhs <- substr(sp_rhs, 1, nchar(sp_rhs) -4)
    sp_rhs <- gsub("G", "G[1:n, 1:n]", sp_rhs)
    sp_rhs <- gsub("G_between", "G_between[1:n, 1:n]", sp_rhs)
    sp_rhs <- gsub("G_within", "G_within[1:n, 1:n]", sp_rhs)
    sp_rhs <- gsub("G_within", "G_within[1:n, 1:n]", sp_rhs)
    sp_rhs <- gsub("G_within_0", "G_within_0[1:n, 1:n]", sp_rhs)
    sp_rhs <- gsub("G_within_1", "G_within_1[1:n, 1:n]", sp_rhs)
    sp_rhs <- gsub("G_within_01", "G_within_01[1:n, 1:n]", sp_rhs)
    sp_rhs <- gsub("G_within_10", "G_within_10[1:n, 1:n]", sp_rhs)
    sp_rhs <- gsub("G_heterogeneity", "G_heterogeneity[1:n, 1:n]", sp_rhs)
    sp_rhs <- gsub("I", "I[1:n, 1:n]", sp_rhs)
    fix_rhs <- paste0(sp_rhs, " %*% ", lin_rhs)
    fix_rhs <- gsub(" \n   ", "", fix_rhs)
  }  else if (model == "model_A") {
    rhs <- rhs.vars(nls_formula_fit)
    lhs <- lhs.vars(nls_formula_fit)
    controls <- 3:length(rhs)
    controls <- controls[controls %% 2 == 0]
    rhs[controls] <- paste0(rhs[controls], "[1:n]")
    sp_rhs <- rhs[2]
    sp_rhs <- gsub("G", "G[1:n, 1:n]", sp_rhs)
    sp_rhs <- gsub("I", "I[1:n, 1:n]", sp_rhs)
    sp_rhs <- gsub("G_heterogeneity", "G_heterogeneity[1:n, 1:n]", sp_rhs)
    sp_rhs <- gsub("Ones", "Ones[1:n]", sp_rhs)

    fix_rhs <- paste0(rhs[1], " * ", sp_rhs, " + ", paste0(rhs[3:length(rhs)], collapse = " + "))
    fix_rhs <- gsub(" \n   ", "", fix_rhs)
  }

  discount_parameter_choice <- function(hypothesis) {
    switch(hypothesis,
           lim = "log(max_egv - (abs(phi) - delta))",
           het = "log(max_egv - (abs(phi) - delta)), log(max_egv - (abs(gamma) - delta))",
           het_l = "log(max_egv - (abs(theta_0) - delta)), log(max_egv - (abs(theta_1) - delta))",
           het_r = "log(max_egv - (abs(eta_0) - delta)), log(max_egv - (abs(eta_1) - delta))",
           par = "log(max_egv - (abs(phi_within) - delta)), log(max_egv - (abs(phi_between) - delta))",
           par_split_with = "log(max_egv - (abs(phi_within_0) - delta)), log(max_egv - (abs(phi_within_1) - delta)), log(max_egv - (abs(phi_between) - delta))",
           par_split_btw = "log(max_egv - (abs(phi_within) - delta)), log(max_egv - (abs(phi_between_01) - delta)), log(max_egv - (abs(phi_between_10) - delta))",
           par_split_with_btw = "log(max_egv - (abs(phi_within_0) - delta)), log(max_egv - (abs(phi_within_1) - delta)), log(max_egv - (abs(phi_between_01) - delta)), log(max_egv - (abs(phi_between_10) - delta))")
  }

  fix_rhs <- paste0("c(", fix_rhs, ",", discount_parameter_choice(hypothesis), ")")
  nls_formula_fit_bou <- paste0(lhs.vars(nls_formula_fit), " ~ ", fix_rhs)
  nls_formula_fit_bou <- formula(nls_formula_fit_bou)

  res <- list(nls_formula_fit = nls_formula_fit_bou, e = e_new, to_weight = new_to_weight)
  return(res)
}
