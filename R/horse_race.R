#' Compare the explanatory power of parameter.dependent network centrality measures with those of standard measures of network centrality.
#'
#' @param formula an object of class \code{formula}: a symbolic description of the model to be fitted. The constant (i.e. intercept) and the autogressive parameter needs not to be specified.
#' @param centralities \ifelse{latex}{at least one of \code{c("indegree","outdegree","degree", "betweenness", "incloseness",} \cr \code{"outcloseness", "closeness", "eigenvector")}.}{at least one of \code{c("indegree","outdegree","degree", "betweenness", "incloseness", "outcloseness", "closeness", "eigenvector")}.}
#' @param directed logical. \code{TRUE} if the social network is directed, \code{FALSE} otherwise.
#' @param weighted logical. \code{TRUE} if the social network is weighted, \code{FALSE} otherwise.
#' @param normalization \ifelse{latex}{Default is NULL. Alternatively, it can be set to \code{c("bygraph","bycomponent",}\cr \code{"bymaxcomponent","bymaxgraph")}. See details.}{Default is NULL. Alternatively, it can be set to \code{c("bygraph","bycomponent","bymaxcomponent","bymaxgraph")}. See details.}
#' @param data an object of class \code{data.frame} containing the variables in the model. If data are longitudinal, observations must be ordered by time period and then by individual.
#' @param unobservables a numeric vector used to obtain an unbiased estimate of the parameter.dependent centrality when the network is endogenous. See details.
#' @param G an object of class \code{Matrix} representing the social network. Row and column names must be specified and match the order of the observations in \code{data}.
#' @param model string. One of \code{c("model_A","model_B")}. See details.
#' @param estimation string. One of \code{c("NLLS","MLE")}. They are used to implement respectively a non-linear least square and a maximum likelihood estimator.
#' @param endogeneity logical. Default is \code{FALSE}. If \code{TRUE}, \code{net_dep} implements a two-step correction procedure to control for the endogeneity of the network.
#' @param first_step \ifelse{latex}{Default is NULL. If \code{endogeneity = TRUE}, it requires to specify one of \code{c("standard",}\cr\code{"fe", "shortest", "coauthors", "degree")}. See details.}{Default is NULL. If \code{endogeneity = TRUE}, it requires to specify one of \code{c("standard","fe", "shortest", "coauthors", "degree")}. See details.}
#' @param data_first_step an optional object of class \code{data.frame}. If provided, it is used to implement the first step of the estimation when \code{endogeneity = TRUE}.
#' @param exclusion_restriction an object of class \code{Matrix} representing the exogenous matrix used to instrument the endogenous social network, if \code{unobservables} is non-\code{NULL}.  Row and column names must be specified and match the order of the observations in \code{data}.
#' @param start.val an optional list containing the starting values for the estimations. Object names must match the names provided in \code{formula}. It is also required to specify the value of both the constant and the decay parameter(s).
#' @param to_weight an optional vector of weights to be used in the fitting process to indicate that different observations have different variances. Should be \code{NULL} or a numeric vector. If non-\code{NULL}, weighted non-linear least squares (if \code{estimation = "NLLS"}) or weighted maximum likelihood  (if \code{estimation = "MLE"}) is estimated.
#' @param time_fixed_effect an optional string. It indicates the name of the time index used in formula. It is used for models with longitudinal data.
#' @param mle_controls a list allowing the user to set upper and lower bounds for control variables in MLE estimation and the variance for the ML estimator.
#' @return A list of two objects:
#' \itemize{
#' \item A list of estimates, each one setting the decay parameter to zero, and adding one of the \code{centralities} to the specification of \code{formula}. The last object adds to \code{formula} all the selected \code{centralities} and the decay parameter is set different from zero.
#' \item An object of class \code{data.frame} containing the computed centrality measures.
#' \item A list of first-step estimations used to correct the effect of centrality measures when the network is endogenous.
#' }
#' @details A number of different normalization are available to the user:
#' \itemize{
#' \item \code{bygraph} and \code{bycomponent} are used to divide \emph{degree} and \emph{closeness} centrality by \eqn{n - 1}, and \emph{betweenness} centrality by \eqn{(n - 1) * (n - 2)} if \code{directed = TRUE}, or by \eqn{(n - 1)*(n - 2)/2} if \code{directed = FALSE}. In the former case (i.e. \code{bygraph}), \emph{n} is equal to the number of nodes in the network In the latter case (i.e. \code{bycomponent}), \emph{n} is equal to the number of nodes of the component in which the node is embedded.
#' \item \code{bymaxgraph} and \code{bymaxcomponent} are used to divide \emph{degree}, \emph{betweenness} and \emph{closeness} centrality by the maximum value of the centrality of the network (\code{bymaxgraph}) or component (\code{bymaxcomponent}) in which the node is embedded.
#' }
#' If the network is endogenous, the user is required to run separately \code{net_dep} and extract from the resulting object the vector of unobservables necessary for obtaining an unbiased estimate of the parameter.dependent centrality. This vector can be passed through the argument \code{unobservables}.\cr
#' If \code{endogeneity = TRUE}, a two-step estimation is implemented to control for network endogeneity. The argument \code{first_step} is used to control for the specification of the first-step model, e.g.:
#' \itemize{
#' \item \code{first_step = "standard"} is used when agents' connection are predicted by the differences in their characteristics (i.e. those on the right hand side of \code{formula}), and an \code{exclusion_restriction}: i.e., their connections in a different network.
#' \item \code{first_step = "fe"} adds individual fixed effects to the \code{standard} model, as in Graham (2017).
#' \item \code{first_step = "shortest"} adds to the \code{standard} model, the shortest distance between \emph{i} and \emph{j}, excluding the link between \emph{i} and \emph{j} itself, as in Fafchamps et al (2010).
#' \item \code{first_step = "coauthor"} adds to the \code{standard} model, the number of shared connections between \emph{i} and \emph{j}, as in Graham (2015).
#' \item \code{first_step = "degree"} adds to the \code{standard} model, the difference in the degree centrality of \emph{i} and \emph{j}.
#' }
#' For additional details, see the vignette.
#' @references
#' Battaglini M., V. Leone Sciabolazza, E. Patacchini, S. Peng (2018), "Econet: An R package for the Estimation of parameter.dependent centrality measures", NBER Working Paper (24442). \cr
#' @seealso \code{\link{net_dep}}
#' @examples
#' \donttest{
#' # Load data
#' data("db_cosponsor")
#' data("G_alumni_111")
#' db_model_B <- subset(db_cosponsor, time  ==  3)
#' G_model_B <- G_cosponsor_111
#' G_exclusion_restriction <- G_alumni_111
#' are_factors <- c("party", "gender", "nchair")
#' db_model_B[are_factors] <- lapply(db_model_B[are_factors], factor)
#'
#' # Specify formula
#' f_model_B <- formula("les ~gender + party + nchair")
#'
#' # Specify starting values
#' starting <- c(alpha = 0.214094,
#'              beta_gender1 = -0.212706,
#'              beta_party1 = 0.478518,
#'              beta_nchair1 = 3.09234,
#'              beta_betweenness = 7.06287e-05,
#'              phi = 0.344787)
#'
#' # Fit model
#' horse_model_B <- horse_race(formula = f_model_B,
#'               centralities = "betweenness",
#'               directed = TRUE, weighted = TRUE,
#'               data = db_model_B, G = G_model_B,
#'               model = "model_B", estimation = "NLLS",
#'               start.val = starting)
#'
#' # Store and print results
#' summary(horse_model_B)
#' summary(horse_model_B, centrality = "betweenness")
#' horse_model_B$centrality
#' }
#' # WARNING, This toy example is provided only for runtime execution.
#' # Please refer to previous examples for sensible calculations.
#' data("db_alumni_test")
#' data("G_model_A_test")
#' db_model <- db_alumni_test
#' G_model <- G_model_A_test
#' f_model <- formula("les ~ dw")
#' horse_model_test <- horse_race(formula = f_model, centralities = "betweenness",
#'                             directed = TRUE, weighted = FALSE, normalization = NULL,
#'                             data = db_model, unobservables = NULL, G = G_model,
#'                             model = "model_A", estimation = "NLLS",
#'                             start.val = c(alpha = -0.31055275,
#'                                           beta_dw = 1.50666982,
#'                                           beta_betweenness = 0.09666742,
#'                                           phi = 16.13035695))
#' summary(horse_model_test)
#' @import spatstat.utils
#' @export
horse_race <- function(
  formula = formula(),
  centralities = c("indegree","outdegree","degree", "betweenness",
                   "incloseness", "outcloseness", "closeness", "eigenvector"),
  directed = FALSE,
  weighted = FALSE,
  normalization = FALSE,
  data = list(),
  unobservables = list(),
  G = list(),
  model = c("model_A","model_B"),
  estimation = c("NLLS","MLE"),
  endogeneity = FALSE,
  first_step = NULL,
  data_first_step = NULL,
  exclusion_restriction = NULL,
  start.val = NULL,
  to_weight = NULL,
  time_fixed_effect = NULL,
  mle_controls = NULL) {

  if (!is.null(time_fixed_effect)) {
    tt <- data[[time_fixed_effect]]
  } else {
    tt <- NULL
  }

  check_1 <- c("indegree","outdegree","degree", "betweenness", "incloseness",
               "outcloseness", "closeness", "eigenvector")
  check_2 <- c("indegree", "outdegree", "incloseness", "outcloseness")
  check_3 <- c("degree", "closeness")

    if (sum(centralities %in% check_1) != length(centralities)) {
      stop('Incorrect field defined')
    }

    if ((sum(centralities %in% check_2) > 0 & (directed == FALSE)) |
       (sum(centralities %in% check_3) > 0 & (directed == TRUE))) {
      stop('Incorrect field defined')
    }

    if (!is.null(exclusion_restriction)) {
      if ((nrow(G) != nrow(exclusion_restriction)) |
         (nrow(G) != nrow(exclusion_restriction))) {
        stop('exclusion restriction must be a squared matrix')
      }
    }

    if (is.null(tt)) {
      if (!is.null(exclusion_restriction)) {
        first_step_cent_mes <- compute_centralities(exclusion_restriction,
                                directed = directed, weighted = weighted,
                                normalization = normalization)
      }
      second_step_cent_mes <- compute_centralities(G,
                                directed = directed, weighted = weighted,
                                normalization = normalization)
    } else {
      k <- unique(tt)
      if (!is.null(exclusion_restriction)) {
      first_step_cent_mes <- Reduce("rbind",lapply(k, function(kk) {
        tmp <- exclusion_restriction[tt == kk, tt == kk]
        compute_centralities(tmp, directed =  directed, weighted = weighted,
                             normalization = normalization)
      }))
      }
      second_step_cent_mes <- Reduce("rbind",lapply(k, function(kk) {
        tmp <- G[tt == kk, tt == kk]
        compute_centralities(tmp, directed =  directed, weighted = weighted,
                             normalization = normalization)
      }))
    }

    if (directed  ==  FALSE) {
      if (!is.null(exclusion_restriction)) {
        colnames(first_step_cent_mes)[colnames(first_step_cent_mes)
                                      == "indegree"] <- "degree"
        colnames(first_step_cent_mes)[colnames(first_step_cent_mes)
                                      == "incloseness"] <- "closeness"
      }
      colnames(second_step_cent_mes)[colnames(second_step_cent_mes)
                                     == "indegree"] <- "degree"
      colnames(second_step_cent_mes)[colnames(second_step_cent_mes)
                                     == "incloseness"] <- "closeness"
    }

    selected_measures <- colnames(second_step_cent_mes)[
                         colnames(second_step_cent_mes) %in% centralities]
    if (!is.null(exclusion_restriction)) {
      first_step_cent_mes <- as.matrix(first_step_cent_mes[,
                              c(selected_measures), drop = FALSE])
      colnames(first_step_cent_mes) <- paste0("instrument_",
                              colnames(first_step_cent_mes))
    }

    second_step_cent_mes <- as.matrix(second_step_cent_mes[,
                              c(selected_measures), drop = FALSE])
    if (!is.null(colnames(second_step_cent_mes))) {
      colnames(second_step_cent_mes) <- centralities
    } else {
      colnames(second_step_cent_mes) <- colnames(second_step_cent_mes)
    }

    data <- cbind(data, second_step_cent_mes)

    if (!is.null(exclusion_restriction)) {
      data <- cbind(data, first_step_cent_mes)
    }

    formula_fit <- formula
    regressors <- paste(termsinformula(formula_fit), collapse = " + ")
    dependent <- as.character(formula_fit)[2]

    if (!is.null(exclusion_restriction)) {
      first_step_list <- list()
      unobs_cent_mes <- list()
      for (i in 1:ncol(first_step_cent_mes)) {
        to_fit <- formula(paste0(colnames(second_step_cent_mes)[i], " ~ ",
                                 regressors, " + ",
                                 colnames(first_step_cent_mes)[i]))
        first_step_list[[i]] <- lm(to_fit, data)
        unobs_cent_mes[[i]] <- residuals(first_step_list[[i]])
      }

      unobs_cent_mes <- as.matrix(Reduce("cbind", unobs_cent_mes))
      colnames(unobs_cent_mes) <- paste0("unobs_", centralities)
      data <- cbind(data, unobs_cent_mes)
    }

    second_step_list <- list()
    for (i in 1:ncol(second_step_cent_mes)) {
      if (!is.null(exclusion_restriction)) {
        to_fit <- formula(paste0(dependent," ~ ", regressors,
                                 " + ", colnames(second_step_cent_mes)[i],
                                 " + ", colnames(unobs_cent_mes)[i]))
      } else {
        to_fit <- formula(paste0(dependent," ~ ", regressors,
                                 " + ", colnames(second_step_cent_mes)[i]))
      }
      second_step_list[[i]] <- lm(to_fit, data, weights = to_weight)
    }

    data_fit <- data
    if (!is.null(exclusion_restriction)) {
      formula_fit <- paste0(dependent," ~ ", regressors, " + ",
                            paste(colnames(second_step_cent_mes),
                                  collapse = " + ")," + ",
                            paste(colnames(unobs_cent_mes),
                                  collapse = " + ")," + unobservables")

      data_fit$unobservables <- unobservables
    } else {
      formula_fit <- paste0(dependent," ~ ", regressors, " + ",
                            paste(colnames(second_step_cent_mes),
                                  collapse = " + "))
    }


    formula_fit <- formula(formula_fit)

    hr <- net_dep(formula = formula_fit, data = data_fit, G = G, model = model,
                  estimation = estimation, hypothesis = "lim",
                  endogeneity = FALSE, correction = "heckman",
                  first_step = NULL, z = NULL, data_first_step = NULL,
                  exclusion_restriction = exclusion_restriction,
                  start.val = start.val, to_weight = to_weight,
                  time_fixed_effect = time_fixed_effect,
                  mle_controls = mle_controls)

    if (!is.null(exclusion_restriction)) {
      second_step_list[[length(second_step_list) + 1]] <- hr[[1]]
      names(second_step_list) <- c(centralities, "parameter.dependent")
      first_step_list[[length(second_step_list) + 1]] <- first_step
      names(first_step_list) <- c(first_step_list, "parameter.dependent")
      second_step_cent_mes <- data.frame(second_step_cent_mes,
                                         "parameter.dependent" = hr[[2]])
      res <- list(second_step = second_step_list,
                  centrality = second_step_cent_mes,
                  first_step = first_step_list)

    }
    else {
      second_step_list[[length(second_step_list) + 1]] <- hr[[1]]
      names(second_step_list) <- c(centralities, "parameter.dependent")
      second_step_cent_mes <- data.frame(second_step_cent_mes,
                                         "parameter.dependent" = hr[[2]])
      res <- list(second_step = second_step_list,
                  centrality = second_step_cent_mes,
                  first_step = NULL)
    }

    class(res) <- "econet"
    attributes(res)$attr <- "horse_race"

    return(res)
  }

