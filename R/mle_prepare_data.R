#' Create a list of objects to be passed to the functions implementing a ML estimator.
#'
#' @param dependent_variable a numeric vector.
#' @param X an object of class \code{data.frame} containing the control variables in the model.
#' @param G G an object of class \code{Matrix} representing the social network.
#' @param z numeric vector. It specifies the source of heterogeneity for peer effects when \code{hypothesis} is equal to \code{"het"}, \code{"het_l"}, or \code{"het_r"}. Alternatively, it specifies the groups in which the network should be partitioned when \code{hypothesis = "par"}.
#' @param e list of objects created either by the user (e.g. \code{G}) or by the function \code{net_dep} (e.g. \emph{"unobservables"})
#' @param start.val a list containing the starting values for the estimations created by the function \code{net_dep}.
#' @param starting.valuesa a list containing the starting values for the estimations provided by the user.
#' @param model string If "model_A", peer effects are modeled as in Battaglini, Patacchini (2018). If "model_B" , peer effects are modeled as in Battaglini, Leone Sciabolazza, Patacchini (2018).
#' @param hypothesis string. It indicates the how peer effects should be modeled.
#' @param endogeneity logical. If \code{TRUE}, it includes the Heckman correction in the \code{data} and the specification, or it instruments the endogenous adjacency matrix.
#' @param correction string. Default is \code{NULL}. If \code{endogeneity = TRUE}, it is required to specify if the main regression should use an instrumental variable ("iv") or Heckman ("heckman") approach.
#' @param tt a vector of positive integers used as time index. It is used for models with longitudinal data.
#' @param mle_controls numeric. It allows the user to set upper and lower bounds for control variables in MLE estimation.
#' @param kappa a normalization level with defult equals 1.
#' @return a list of two objects: i) a list of starting values; ii) a list of objects containing the data to be used in the regression.
#' @keywords internal
#' @noRd
mle_prepare_data<- function(X, G, z, e, start.val, starting.values, model,
                            hypothesis, endogeneity, correction, tt,
                            mle_controls, kappa) {

  if (endogeneity == TRUE) {
    if (correction == "heckman") {
      X <-data.frame(X, unobservables = e[["unobservables"]])
      if (is.null(start.val)) {
        starting.values[["beta_unobservables"]] <- 0.01
      }
    }
  }

  if (model == "model_A") {

    if (hypothesis == "lim") {

      if (is.null(start.val)) {
        tmp1 <- starting.values
        tmp1[["phi"]] <- 0.01
        tmp1[["alpha"]] <- 0.01
        tmp1[["sigma"]] <- 0.1
        tmp1 <- unlist(tmp1)
      } else {
        tmp1 <- unlist(start.val)
      }

      if (is.null(mle_controls)) {
        tmp2 <- tmp1 + 100
        tmp3 <- tmp1 - 100
        tmp2[names(tmp1) %in% "sigma"] <- 10
        tmp3[names(tmp1) %in% "sigma"] <- 0.01
      } else {
        tmp2 <- tmp1 + mle_controls[[1]]
        tmp3 <- tmp1 - mle_controls[[1]]
        tmp2[names(tmp1) %in% "sigma"] <- mle_controls[[2]][1]
        tmp3[names(tmp1) %in% "sigma"] <- mle_controls[[2]][2]
      }

      tmp2[names(tmp1) %in% "phi"] <- 1 / kappa
      tmp3[names(tmp1) %in% "phi"] <- -1 / kappa

    e[["G"]] <- G

    } else if (hypothesis == "het") {

      if (is.null(start.val)) {
        tmp1 <- starting.values
        tmp1[["alpha"]] <- 0.01
        tmp1[["phi"]] <- 0.01
        tmp1[["gamma"]] <- 0.01
        tmp1[["sigma"]] <- 0.1
        tmp1 <- unlist(tmp1)
      } else {
        tmp1 <- unlist(start.val)
      }

      if (is.null(mle_controls)) {
        tmp2 <- tmp1 + 100
        tmp3 <- tmp1 - 100
        tmp2[names(tmp1) %in% "sigma"] <- 10
        tmp3[names(tmp1) %in% "sigma"] <- 0.01
      } else {
        tmp2 <- tmp1 + mle_controls[[1]]
        tmp3 <- tmp1 - mle_controls[[1]]
        tmp2[names(tmp1) %in% "sigma"] <- mle_controls[[2]][1]
        tmp3[names(tmp1) %in% "sigma"] <- mle_controls[[2]][2]
      }

      tmp2[names(tmp1) %in% "phi"] <- 1 / kappa
      tmp3[names(tmp1) %in% "phi"] <- -1 / kappa

      e[["G_heterogeneity"]] <- diag(z)
    }

  } else if (model == "model_B") {

    if (hypothesis == "lim") {
      if (is.null(start.val)) {
        tmp1 <- starting.values
        tmp1[["phi"]] <- 0.01
        tmp1[["alpha"]] <- 0.01
        tmp1[["sigma"]] <- 0.1
        tmp1 <- unlist(tmp1)
      } else {
        tmp1 <- unlist(start.val)
      }

      if (is.null(mle_controls)) {
        tmp2 <- tmp1 + 100
        tmp3 <- tmp1 - 100
        tmp2[names(tmp1) %in% "sigma"] <- 10
        tmp3[names(tmp1) %in% "sigma"] <- 0.01
      } else {
        tmp2 <- tmp1 + mle_controls[[1]]
        tmp3 <- tmp1 - mle_controls[[1]]
        tmp2[names(tmp1) %in% "sigma"] <- mle_controls[[2]][1]
        tmp3[names(tmp1) %in% "sigma"] <- mle_controls[[2]][2]
      }

      tmp2[names(tmp1) %in% "phi"] <- 1 / kappa
      tmp3[names(tmp1) %in% "phi"] <- -1 / kappa

      e[["G"]] <- G

    } else if (hypothesis == "het_l") {

      e[["G"]] <- G
      e[["G_heterogeneity"]] <- diag(z)

      if (is.null(start.val)) {
        tmp1 <- starting.values
        tmp1[["alpha"]] <- 0.01
        tmp1[["theta_0"]] <- 0.01
        tmp1[["theta_1"]] <- 0.01
        tmp1[["sigma"]] <- 0.1
        tmp1 <- unlist(tmp1)
      } else {
        tmp1 <- unlist(start.val)

      }

      if (is.null(mle_controls)) {
        tmp2 <- tmp1 + 100
        tmp3 <- tmp1 - 100
        tmp2[names(tmp1) %in% "sigma"] <- 10
        tmp3[names(tmp1) %in% "sigma"] <- 0.01
      } else {
        tmp2 <- tmp1 + mle_controls[[1]]
        tmp3 <- tmp1 - mle_controls[[1]]
        tmp2[names(tmp1) %in% "sigma"] <- mle_controls[[2]][1]
        tmp3[names(tmp1) %in% "sigma"] <- mle_controls[[2]][2]
      }

      tmp2[names(tmp1) %in% "theta_1"] <- 2 / kappa
      tmp3[names(tmp1) %in% "theta_1"] <- -2 / kappa

      tmp2[names(tmp1) %in% "theta_0"] <- 1 / kappa
      tmp3[names(tmp1) %in% "theta_0"] <- -1 / kappa

    } else if (hypothesis == "het_r") {

      e[["G"]] <- G
      e[["G_heterogeneity"]] <- diag(z)

      if (is.null(start.val)) {
        tmp1 <- starting.values
        tmp1[["alpha"]] <- 0.01
        tmp1[["eta_0"]] <- 0.01
        tmp1[["eta_1"]] <- 0.01
        tmp1[["sigma"]] <- 0.1
        tmp1 <- unlist(tmp1)
      } else {
        tmp1 <- unlist(start.val)
      }

      if (is.null(mle_controls)) {
        tmp2 <- tmp1 + 100
        tmp3 <- tmp1 - 100
        tmp2[names(tmp1) %in% "sigma"] <- 10
        tmp3[names(tmp1) %in% "sigma"] <- 0.01
      } else {
        tmp2 <- tmp1 + mle_controls[[1]]
        tmp3 <- tmp1 - mle_controls[[1]]
        tmp2[names(tmp1) %in% "sigma"] <- mle_controls[[2]][1]
        tmp3[names(tmp1) %in% "sigma"] <- mle_controls[[2]][2]
      }

      tmp2[names(tmp1) %in% "theta_1"] <- 2 / kappa
      tmp3[names(tmp1) %in% "theta_1"] <- -2 / kappa

      tmp2[names(tmp1) %in% "theta_0"] <- 1 / kappa
      tmp3[names(tmp1) %in% "theta_0"] <- -1 / kappa

    } else if (hypothesis == "par") {

      if (is.null(tt)) {
        partitions <- party_connection(z, G)
      } else {
        partitions <- time_party_connection(z, G, tt)
      }

      e[["G_within"]] <- partitions[[1]]
      e[["G_between"]] <- partitions[[2]]

      if (is.null(start.val)) {
        tmp1 <- starting.values
        tmp1[["alpha"]] <- 0.01
        tmp1[["phi_within"]] <- 0.01
        tmp1[["phi_between"]] <- 0.01
        tmp1[["sigma"]] <- 0.1
        tmp1 <- unlist(tmp1)
      } else {
        tmp1 <- unlist(start.val)
      }

      if (is.null(mle_controls)) {
        tmp2 <- tmp1 + 100
        tmp3 <- tmp1 - 100
        tmp2[names(tmp1) %in% "sigma"] <- 10
        tmp3[names(tmp1) %in% "sigma"] <- 0.01
      } else {
        tmp2 <- tmp1 + mle_controls[[1]]
        tmp3 <- tmp1 - mle_controls[[1]]
        tmp2[names(tmp1) %in% "sigma"] <- mle_controls[[2]][1]
        tmp3[names(tmp1) %in% "sigma"] <- mle_controls[[2]][2]
      }
      tmp2[names(tmp1) %in% "phi_within"] <- 1 / max(e[["G_within"]])
      tmp3[names(tmp1) %in% "phi_within"] <- -1 / max(e[["G_within"]])
      tmp2[names(tmp1) %in% "phi_between"] <- 1 / max(e[["G_between"]])
      tmp3[names(tmp1) %in% "phi_between"] <- -1 / max(e[["G_between"]])
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
        tmp1 <- starting.values
        tmp1[["alpha"]] <- 0.01
        tmp1[["phi_within_0"]] <- 0.01
        tmp1[["phi_within_1"]] <- 0.01
        tmp1[["phi_between"]] <- 0.01
        tmp1[["sigma"]] <- 0.1
        tmp1 <- unlist(tmp1)
      } else {
        tmp1 <- unlist(start.val)
      }

      if (is.null(mle_controls)) {
        tmp2 <- tmp1 + 100
        tmp3 <- tmp1 - 100
        tmp2[names(tmp1) %in% "sigma"] <- 10
        tmp3[names(tmp1) %in% "sigma"] <- 0.01
      } else {
        tmp2 <- tmp1 + mle_controls[[1]]
        tmp3 <- tmp1 - mle_controls[[1]]
        tmp2[names(tmp1) %in% "sigma"] <- mle_controls[[2]][1]
        tmp3[names(tmp1) %in% "sigma"] <- mle_controls[[2]][2]
      }
      tmp2[names(tmp1) %in% "phi_within_0"] <- 1 / max(e[["G_within_0"]])
      tmp3[names(tmp1) %in% "phi_within_0"] <- -1 / max(e[["G_within_0"]])
      tmp2[names(tmp1) %in% "phi_within_1"] <- 1 / max(e[["G_within_1"]])
      tmp3[names(tmp1) %in% "phi_within_1"] <- -1 / max(e[["G_within_1"]])
      tmp2[names(tmp1) %in% "phi_between"] <- 1 / max(e[["G_between"]])
      tmp3[names(tmp1) %in% "phi_between"] <- -1 / max(e[["G_between"]])
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
        tmp1 <- starting.values
        tmp1[["alpha"]] <- 0.01
        tmp1[["phi_within"]] <- 0.01
        tmp1[["phi_between_01"]] <- 0.01
        tmp1[["phi_between_10"]] <- 0.01
        tmp1[["sigma"]] <- 0.1
        tmp1 <- unlist(tmp1)
      } else {
        tmp1 <- unlist(start.val)
      }

      if (is.null(mle_controls)) {
        tmp2 <- tmp1 + 100
        tmp3 <- tmp1 - 100
        tmp2[names(tmp1) %in% "sigma"] <- 10
        tmp3[names(tmp1) %in% "sigma"] <- 0.01
      } else {
        tmp2 <- tmp1 + mle_controls[[1]]
        tmp3 <- tmp1 - mle_controls[[1]]
        tmp2[names(tmp1) %in% "sigma"] <- mle_controls[[2]][1]
        tmp3[names(tmp1) %in% "sigma"] <- mle_controls[[2]][2]
      }
      tmp2[names(tmp1) %in% "phi_within"] <- 1 / max(e[["G_within"]])
      tmp3[names(tmp1) %in% "phi_within"] <- -1 / max(e[["G_within"]])
      tmp2[names(tmp1) %in% "phi_between_01"] <- 1 / max(e[["G_between_01"]])
      tmp3[names(tmp1) %in% "phi_between_01"] <- -1 / max(e[["G_between_01"]])
      tmp2[names(tmp1) %in% "phi_between_10"] <- 1 / max(e[["G_between_10"]])
      tmp3[names(tmp1) %in% "phi_between_10"] <- -1 / max(e[["G_between_10"]])
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
        tmp1 <- starting.values
        tmp1[["alpha"]] <- 0.01
        tmp1[["phi_within_0"]] <- 0.01
        tmp1[["phi_within_1"]] <- 0.01
        tmp1[["phi_between_01"]] <- 0.01
        tmp1[["phi_between_10"]] <- 0.01
        tmp1[["sigma"]] <- 0.1
        tmp1 <- unlist(tmp1)
      } else {
        tmp1 <- unlist(start.val)
      }

      if (is.null(mle_controls)) {
        tmp2 <- tmp1 + 100
        tmp3 <- tmp1 - 100
        tmp2[names(tmp1) %in% "sigma"] <- 10
        tmp3[names(tmp1) %in% "sigma"] <- 0.01
      } else {
        tmp2 <- tmp1 + mle_controls[[1]]
        tmp3 <- tmp1 - mle_controls[[1]]
        tmp2[names(tmp1) %in% "sigma"] <- mle_controls[[2]][1]
        tmp3[names(tmp1) %in% "sigma"] <- mle_controls[[2]][2]
      }
      tmp2[names(tmp1) %in% "phi_within_0"] <- 1 / max(e[["G_within_0"]])
      tmp3[names(tmp1) %in% "phi_within_0"] <- -1 / max(e[["G_within_0"]])
      tmp2[names(tmp1) %in% "phi_within_1"] <- 1 / max(e[["G_within_1"]])
      tmp3[names(tmp1) %in% "phi_within_1"] <- -1 / max(e[["G_within_1"]])
      tmp2[names(tmp1) %in% "phi_between_01"] <- 1 / max(e[["G_between_01"]])
      tmp3[names(tmp1) %in% "phi_between_01"] <- -1 / max(e[["G_between_01"]])
      tmp2[names(tmp1) %in% "phi_between_10"] <- 1 / max(e[["G_between_10"]])
      tmp3[names(tmp1) %in% "phi_between_10"] <- -1 / max(e[["G_between_10"]])
    }
  }

  starting.values <- tmp1
  boundU <- tmp2
  boundL <- tmp3

  res <- list(starting.values, boundU, boundL, e)

  return(res)

}
