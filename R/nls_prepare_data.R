#' Create a list of objects to be passed to the function \code{nlsLM}
#'
#' @param dependent_variable a numeric vector.
#' @param X an object of class \code{data.frame} containing the control variables in the model.
#' @param G G an object of class \code{Matrix} representing the social network.
#' @param z numeric vector. It specifies the source of heterogeneity for peer effects when \code{hypothesis} is equal to \code{"het"}, \code{"het_l"}, or \code{"het_r"}. Alternatively, it specifies the groups in which the network should be partitioned when \code{hypothesis = "par"}.
#' @param e list of objects created either by the user (e.g. \code{G}) or by the function \code{net_dep} (e.g. \emph{"unobservables"})
#' @param start.val a list containing the starting values for the estimations created by the function \code{net_dep}.
#' @param starting.values a list containing the starting values for the estimations provided by the user.
#' @param model string If "model_A", peer effects are modeled as in Battaglini, Patacchini (2018). If "model_B" , peer effects are modeled as in Battaglini, Leone Sciabolazza, Patacchini (2018).
#' @param hypothesis string. It indicates the how peer effects should be modeled.
#' @param endogeneity logical. If \code{TRUE}, it includes the Heckman correction in the \code{data} and the specification, or it instruments the endogenous adjacency matrix.
#' @param correction string. Default is \code{NULL}. If \code{endogeneity = TRUE}, it is required to specify if the main regression should use an instrumental variable ("iv") or Heckman ("heckman") approach.
#' @param tt a vector of positive integers used as time index. It is used for models with longitudinal data.
#' @param nls_controls a list allowing the user to set upper and lower bounds for control variables in NLLS estimation.
#' @return a list of three objects: i) an object of class \code{formula}; ii) a list of starting values; iii) a list of objects containing the data to be used in the regression.
#' @keywords internal
#' @noRd
nls_prepare_data <- function(dependent_variable, X, G, z, e, start.val, starting.values, model, hypothesis, endogeneity, correction, tt, nls_controls) {

  parameter_choice <- function(hypothesis) {
    switch(hypothesis,
           lim = "solve_block(I - phi * G)",
           het = "solve_block(I - G %*% (phi * I + gamma * G_heterogeneity))",
           het_l = "solve_block(I - (theta_0 * I - theta_1 * G_heterogeneity) %*% G)",
           het_r = "solve_block(I - G %*% (eta_0 * I - eta_1 * G_heterogeneity))",
           par = "solve_block(I - phi_within * G_within - phi_between * G_between)",
           par_split_with = "solve_block(I - phi_within_0 * G_within_0 - phi_within_1 * G_within_1 - phi_between * G_between)",
           par_split_btw = "solve_block(I - phi_within * G_within - phi_between_01 * G_between_01 - phi_between_10 * G_between_10)",
           par_split_with_btw = "solve_block(I - phi_within_0 * G_within_0 - phi_within_1 * G_within_1 - phi_between_01 * G_between_01 - phi_between_10 * G_between_10)")
  }

  model_choice <- function(x, y, z, model) {
    switch(model,
           model_A = formula(paste0(x, " ~ alpha *", y, "%*% Ones +", z)),
           model_B = formula(paste0(x, " ~ ", y, "%*% (alpha * Ones +", z, ")")))}

  regressors <- paste(paste0("beta_", colnames(X)[-1])," * ", colnames(X)[-1],
                      collapse = " + ")

  sel <- which(names(starting.values) %in% colnames(X))
  names(starting.values)[sel] <- paste0("beta_", names(starting.values)[sel])

  if (endogeneity == TRUE) {
    if (correction == "heckman") {
      formula_regressors <- paste0(regressors,
                                   " + beta_unobservables * unobservables")
      if (is.null(start.val)) {
        starting.values[["beta_unobservables"]] <- 0.01
      } else {
        names(starting.values)[names(starting.values) %in%
                                 "unobservables"] <- "beta_unobservables"
      }
    } else {
      formula_regressors <- regressors
    }
  } else {
    formula_regressors <- regressors
  }

  formula_parameter_dependent <- parameter_choice(hypothesis)

  formula_fit <- model_choice(x = dependent_variable,
                              y = formula_parameter_dependent,
                              z= formula_regressors, model)

  if (hypothesis == "lim") {
    e[["G"]] <- G

    if (is.null(start.val)) {
      starting.values[["phi"]] <- 0.01
    }
  } else if (hypothesis == "het") {
    e[["G"]] <- G
    e[["G_heterogeneity"]] <- diag(z)

    if (is.null(start.val)) {
      starting.values[["phi"]] <- 0.01
      starting.values[["gamma"]] <- 0.01
    }
  } else if (hypothesis == "het_l" | hypothesis == "het_r") {

    e[["G"]] <- G
    e[["G_heterogeneity"]] <- diag(z)

    if (is.null(start.val) & hypothesis == "het_l") {
      starting.values[["theta_0"]] <- 0.01
      starting.values[["theta_1"]] <- 0.01
    } else if (is.null(start.val) & hypothesis == "het_r") {
      starting.values[["eta_0"]] <- 0.01
      starting.values[["eta_1"]] <- 0.01
    }
  } else if (hypothesis == "par") {

    if (is.null(tt)) {
      partitions <- party_connection(z, G)
    } else {
      partitions <- time_party_connection(z, G, tt)
    }

    e[["G_within"]] <- partitions[[1]]
    e[["G_between"]] <- partitions[[2]]

    if (is.null(start.val)) {
      starting.values[["phi_within"]] <- 0.01
      starting.values[["phi_between"]] <- 0.01
    }
  } else if (hypothesis == "par_split_with") {

    if (is.null(tt)) {
      partitions <- party_connection_split(z, G)
    } else {
      partitions <- time_party_connection_split(z, G, tt)
    }

    e[["G_within_0"]] <- partitions[[1]]
    e[["G_within_1"]] <- partitions[[2]]
    e[["G_between"]] <- partitions[[3]] + partitions[[4]]

    if (is.null(start.val)) {
      starting.values[["phi_within_0"]] <- 0.01
      starting.values[["phi_within_1"]] <- 0.01
      starting.values[["phi_between"]] <- 0.01
    }

  } else if (hypothesis == "par_split_btw") {

    if (is.null(tt)) {
      partitions <- party_connection_split(z, G)
    } else {
      partitions <- time_party_connection_split(z, G, tt)
    }

    e[["G_within"]] <- partitions[[1]] + partitions[[2]]
    e[["G_between_01"]] <- partitions[[3]]
    e[["G_between_10"]] <- partitions[[4]]

    if (is.null(start.val)) {
      starting.values[["phi_within"]] <- 0.01
      starting.values[["phi_between_01"]] <- 0.01
      starting.values[["phi_between_10"]] <- 0.01
    }

  } else if (hypothesis == "par_split_with_btw") {

    if (is.null(tt)) {
      partitions <- party_connection_split(z, G)
    } else {
      partitions <- time_party_connection_split(z, G, tt)
    }

    e[["G_within_0"]] <- partitions[[1]]
    e[["G_within_1"]] <- partitions[[2]]
    e[["G_between_01"]] <- partitions[[3]]
    e[["G_between_10"]] <- partitions[[4]]

    if (is.null(start.val)) {
      starting.values[["phi_within_0"]] <- 0.01
      starting.values[["phi_within_1"]] <- 0.01
      starting.values[["phi_between_01"]] <- 0.01
      starting.values[["phi_between_10"]] <- 0.01
    }
  }

  list(formula_fit, starting.values, e)
}

