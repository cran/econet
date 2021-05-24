#' Implement a number of modifications to the linear-in-means model to obtain different weighted versions of Katz-Bonacich centrality.
#'
#' @param formula an object of class \code{formula}: a symbolic description of the model to be fitted. The constant (i.e. intercept) and the autogressive parameter needs not to be specified.
#' @param data an object of class \code{data.frame} containing the variables in the model. If data are longitudinal, observations must be ordered by time period and then by individual.
#' @param G an object of class \code{Matrix} representing the social network. Row and column names must be specified and match the order of the observations in \code{data}.
#' @param model string. One of \code{c("model_A","model_B")}. See details.
#' @param estimation string. One of \code{c("NLLS","MLE")}. They are used to implement respectively a non-linear least square and a Maximum Likelihood estimator.
#' @param hypothesis string. One of \code{c("lim","het", "het_l", "het_r", "par", "par_split_with", "par_split_btw", "par_split_with_btw")}. See details.
#' @param endogeneity logical. Default is \code{FALSE}. If \code{TRUE}, \code{net_dep} implements a two-step correction procedure to control for the endogeneity of the network.
#' @param correction Default is \code{NULL}. If \code{endogeneity = TRUE}, it is required to specify if the main regression should use an instrumental variable ("iv") or Heckman ("heckman") approach.
#' @param first_step \ifelse{latex}{Default is NULL. If \code{endogeneity = TRUE}, it requires to specify one of \code{c("standard"}\cr\code{"fe", "shortest", "coauthors", "degree")}. See details.}{Default is NULL. If \code{endogeneity = TRUE}, it requires to specify one of \code{c("standard","fe", "shortest", "coauthors", "degree")}. See details.}
#' @param z numeric vector. It specifies the source of heterogeneity for peer effects when \code{hypothesis} is equal to \code{"het"}, \code{"het_l"}, or \code{"het_r"}. Alternatively, it specifies the groups in which the network should be partitioned when \code{hypothesis} is equal to \code{"par"}, \code{"par_split_with"}, \code{"par_split_btw"}, or \code{"par_split_with_btw")}. See details.
#' @param formula_first_step an optional object of class \code{formula}. If provided, it is used to implement the first step of the estimation when \code{endogeneity = TRUE}. The name of the dependent variable must be the same used in the field  \code{formula}.
#' @param exclusion_restriction an object of class \code{Matrix} representing the exogenous matrix used to instrument the endogenous social network, if \code{endogeneity = TRUE}.  Row and column names must be specified and match the order of the observations in \code{data}.
#' @param start.val an optional list containing the starting values for the estimations. Object names must match the names provided in \code{formula}. It is also required to specify the value of both the constant and the decay parameter(s). See details.
#' @param to_weight an optional vector of weights to be used in the fitting process to indicate that different observations have different variances. Should be \code{NULL} or a numeric vector. If non-\code{NULL}, it can be used to fit a weighted non-linear least squares (\code{estimation = "NLLS"}).
#' @param time_fixed_effect an optional string. It indicates the name of the time index used in formula. It is used for models with longitudinal data.
#' @param ind_fixed_effect an optional string. Default is \code{NULL}. It indicates the name of the individual index contained in the data. If provided, individual fixed effects are automatically added to the \code{formula} of the main equation. If \code{endogeneity = TRUE}, the field \code{first_step} is overridden, and automatically set equal to \code{"fe"}. It is used for models with longitudinal data. Observe that inclusion of individual fixed effects is in its beta version. When \code{estimation == "MLE"}, \code{net_dep} is guaranteed to work only if \code{hypothesis == "lim"}.
#' @param mle_controls a list allowing the user to set upper and lower bounds for control variables in MLE estimation and the variance for the ML estimator. See details.
#' @param kappa a normalization level with default equals 1 used in MLE estimation.
#' @param delta Default is \code{NULL}. To be used when \code{estimation = "NLLS"}. It has to be a number between zero (included) and one (excluded). When used, \code{econet} performs a constrained NLLS estimation. In this case, the estimated peer effect parameter, taken in absolute value, is forced to be between the spectral radius of \code{G} and its opposite value. Specifically, \code{delta} is a penalizing factor, decreasing the goodness of fit of the NLLS estimation, when the peer effect parameter approaches one of the two bounds. Observe that very high values of \code{delta} may cause NLLS estimation not to converge.
#' @return A list of three objects: i) Estimates of the main regression; ii) The vector of agents' parameter-dependent centrality; iii) Estimates of the first-step regression (if \code{endogeneity = TRUE})
#' @details Agent's parameter-dependent centrality is obtained as a function of \itemize{
#' \item the agent's characteristics and the performance of its socially connected peers, as in Battaglini, Leone Sciabolazza, Patacchini (2020), if \code{model = "model_B"};
#' \item \ifelse{latex}{the performance of its socially connected peers, as in Battaglini, Patacchini (2018), if \cr \code{model = "model_A"}.}{the performance of its socially connected peers, as in Battaglini, Patacchini (2018), if \code{model = "model_A"}.}
#' }
#' Peer effects are assumed to be homogenous if \code{hypothesis = "lim"}. They are assumed to be heterogenous by setting:\itemize{
#' \item \ifelse{latex}{\code{hypothesis = "het"}, when peers' performance is susceptible to agent's characteristics and \cr \code{model = "model_A"}.}{\code{hypothesis = "het"}, when peers' performance is susceptible to agent's characteristics and \code{model = "model_A"}.}
#' \item \code{hypothesis = "het_l"}, when peers' performance is susceptible to agent's characteristics and \code{model = "model_B"}.
#' \item \code{hypothesis = "het_r"}, when agent's performance is susceptible to peers' characteristics and \code{model = "model_B"}.
#' \item \code{hypothesis = "par"}, when \code{model = "model_B"}, if the network is formed by interactions between and within two different groups.
#' \item \code{hypothesis = "par_split_with"}, when \code{model = "model_B"}, if the network is formed by interactions between and within two different groups, and interactions within each group are different from the other.
#' \item \code{hypothesis = "par_split_btw"}, when \code{model = "model_B"}, if the network is formed by interactions between and within two different groups, and interactions between groups are different according to their direction.
#' \item \code{hypothesis = "par_split_with_btw"}, when \code{model = "model_B"}, if the network is formed by interactions between and within two different groups, interactions within each group are different from the other, and interactions between groups are different according to their direction.
#' }
#' When \code{hypothesis} is equal to \code{"het"}, \code{"het_l"}, or \code{"het_r"}, the argument \code{z} is used to specify the source of heterogeneity: i.e. the attribute affecting the ability of the agent to influence or be influenced by peers.
#' When \code{hypothesis} is equal to \code{"par"}, \code{"par_split_with"}, \code{"par_split_btw"}, or \code{"par_split_with_btw"} the argument \code{"z"} is used to partition observations in two groups: e.g. the generic element \emph{i} of vector \code{z} takes the value 1 if agent \emph{i} is member of the first group, and it takes 2 otherwise. \cr
#' \cr
#' If \code{endogeneity = TRUE}, a two-step estimation is implemented to control for network endogeneity. The argument \code{first_step} is used to control for the specification of the first-step model, e.g.:
#' \itemize{
#' \item \code{first_step = "standard"} is used when agents' connection are predicted by the differences in their characteristics (i.e. those on the right hand side of \code{formula}), and an \code{exclusion_restriction}: i.e., their connections in a different network.
#' \item \code{first_step = "fe"} adds to the \code{standard} model, individual fixed effects, as in Graham (2017).
#' \item \code{first_step = "shortest"} adds to the \code{standard} model, the shortest distance between \emph{i} and \emph{j}, excluding the link between \emph{i} and \emph{j} itself, as in Fafchamps et al (2010).
#' \item \code{first_step = "coauthor"} adds to the \code{standard} model, the number of shared connections between \emph{i} and \emph{j}, as in Graham (2015).
#' \item \code{first_step = "degree"} adds to the \code{standard} model, the difference in the degree centrality of \emph{i} and \emph{j}.
#' }
#' The argument \code{start.val} is used to specify starting estimates. This can be done with a named list. If a factor is present, a value for each \code{treatment contrast} must be provided. Labels of \code{treatment contrast}s must be assigned following R model design standards: e.g., a number is appended to contrast names as in \code{contrasts()}. \cr
#' The starting value referring to the intercept (constant) must be labelled as "alpha". The label(s) for decay parameter(s) must be:
#' \itemize{
#' \item \code{"phi"}, if \code{hypothesis="lim"} or  \code{hypothesis="het"}
#' \item \code{"theta_0","theta_1"}, if \code{hypothesis="het_l"}
#' \item \code{"eta_0","eta_1"}, if \code{hypothesis="het_r"}
#' \item \code{"phi_within","phi_between"}, if \code{hypothesis="par"}
#' \item \code{"phi_within_0","phi_within_1","phi_between"}, if \code{hypothesis="par_split_with"}
#' \item \code{"phi_within","phi_between_0","phi_between_1"}, if \code{hypothesis="par_split_btw"}
#' \item \code{"phi_within_0","phi_within_1","phi_between_0","phi_between_1"}, if \code{hypothesis="par_split_with_btw"}
#' }
#' The interaction term when \code{hypothesis="het"} must be labelled  \code{"gamma"}. The label to be used for unobservables when \code{endogeneity = TRUE} is \code{"unobservables"}. When \code{estimation = "MLE"}, it is required to set also the starting value for the variance of the ML estimator. This should be labelled as "sigma". \cr
#' The argument \code{mle_controls} takes a list of two objects. The first is a named numeric vector used to set upper and lower bounds for control variables. The second object is a vector used to set upper (first value) and lower (second value) bounds for the variance of the Maximum Likelihood estimator.\cr
#' Names in \code{mle_controls} must be equal to those used in \code{start.val}.
#' For additional details, see the vignette.
#' @references
#' Battaglini M., E. Patacchini (2018), "Influencing Connected Legislators," Journal of Political Economy, 126(6): 2277-2322. \cr
#' Battaglini M., V. Leone Sciabolazza, E. Patacchini (2020), "Effectiveness of Connected Legislators,"  American Journal of Political Science, forthcoming. \cr
#' Battaglini M., V. Leone Sciabolazza, E. Patacchini, S. Peng (2020), "Econet: An R package for the Estimation of parameter-dependent centrality measures", Mimeo. \cr
#' Fafchamps, M., M. J. Leij and S. Goyal (2010), “Matching and network effects,” Journal of the European Economic Association, 8(1): 203-231. \cr
#' Graham B. (2015), “Methods of identification in social networks,” Annual Review of Economics, 7, 465 - 485. \cr
#' Graham B. (2017), “An econometric model of network formation with degree heterogeneity,” Econometrica 85 (4), 1033 - 1063. \cr
#' @examples
#' \donttest{
#' # Model A
#'
#' # Load data
#' data("a_db_alumni")
#' data("a_G_alumni_111")
#' db_model_A <- a_db_alumni
#' G_model_A <- a_G_alumni_111
#' are_factors <- c("party", "gender", "nchair", "isolate")
#' db_model_A[are_factors] <- lapply(db_model_A[are_factors] ,factor)
#' db_model_A$PAC <- db_model_A$PAC/1e+06
#'
#' # Specify formula
#' f_model_A <- formula("PAC ~ gender + party + nchair + isolate")
#'
#' # Specify starting values
#' starting <- c(alpha = 0.47325,
#'               beta_gender1 = -0.26991,
#'               beta_party1 = 0.55883,
#'               beta_nchair1 = -0.17409,
#'               beta_isolate1 = 0.18813,
#'               phi = 0.21440)
#'
#' # Fit Linear-in-means model
#' lim_model_A <- net_dep(formula = f_model_A, data = db_model_A,
#'                        G = G_model_A, model = "model_A", estimation = "NLLS",
#'                        hypothesis = "lim", start.val = starting)
#'
#' summary(lim_model_A)
#' lim_model_A$centrality
#'
#' # Test Heterogeneity
#'
#' # Heterogeneous factor
#' z <- as.numeric(as.character(db_model_A$gender))
#'
#' # Specify formula
#' f_het_model_A <- formula("PAC ~ party + nchair + isolate")
#'
#' # Specify starting values
#' starting <- c(alpha = 0.44835,
#'               beta_party1 = 0.56004,
#'               beta_nchair1 = -0.16349,
#'               beta_isolate1 = 0.21011,
#'               beta_z = -0.26015,
#'               phi = 0.34212,
#'               gamma = -0.49960)
#'
#' # Fit model
#' het_model_A <- net_dep(formula = f_het_model_A, data = db_model_A,
#'                        G = G_model_A, model = "model_A", estimation = "NLLS",
#'                        hypothesis = "het", z = z, start.val = starting)
#'
#' summary(het_model_A)
#' het_model_A$centrality
#'
#' # Model B
#'
#' # Load data
#' data("db_cosponsor")
#' data("G_alumni_111")
#' db_model_B <- db_cosponsor
#' G_model_B <- G_cosponsor_111
#' G_exclusion_restriction <- G_alumni_111
#' are_factors <- c("party", "gender", "nchair")
#' db_model_B[are_factors] <- lapply(db_model_B[are_factors], factor)
#'
#' # Specify formula
#' f_model_B <- formula("les ~ gender + party + nchair")
#'
#' # Specify starting values
#' starting <- c(alpha = 0.23952,
#'               beta_gender1 = -0.22024,
#'               beta_party1 = 0.42947,
#'               beta_nchair1 = 3.09615,
#'               phi = 0.40038,
#'               unobservables = 0.07714)
#'
#' # Fit Linear-in-means model
#' lim_model_B <- net_dep(formula = f_model_B, data = db_model_B,
#'                        G = G_model_B, model = "model_B", estimation = "NLLS",
#'                        hypothesis = "lim", endogeneity = TRUE,
#'                        correction = "heckman", first_step = "standard",
#'                        exclusion_restriction = G_exclusion_restriction,
#'                        start.val = starting)
#'
#' summary(lim_model_B)
#' lim_model_B$centrality
#' summary(lim_model_B, print = "first.step")
#'
#' # Test Heterogeneity
#'
#' # Heterogeneous factor (node -level)
#' z <- as.numeric(as.character(db_model_B$gender))
#'
#' # Specify formula
#' f_het_model_B <- formula("les ~ party + nchair")
#'
#' # Specify starting values
#' starting <- c(alpha = 0.23952,
#'               beta_party1 = 0.42947,
#'               beta_nchair1 = 3.09615,
#'               beta_z = -0.12749,
#'               theta_0 = 0.42588,
#'               theta_1 = 0.08007)
#'
#' # Fit model
#' het_model_B_l <- net_dep(formula = f_het_model_B,
#'                          data = db_model_B,
#'                          G = G_model_B, model = "model_B", estimation = "NLLS",
#'                          hypothesis = "het_l", z = z, start.val = starting)
#'
#' # Store and print results
#' summary(het_model_B_l)
#' het_model_B_l$centrality
#'
#' # Specify starting values
#' starting <- c(alpha = 0.04717,
#'               beta_party1 = 0.51713,
#'               beta_nchair1 = 3.12683,
#'               beta_z = 0.01975,
#'               eta_0 = 1.02789,
#'               eta_1 = 2.71825)
#'
#' # Fit model
#' het_model_B_r <- net_dep(formula = f_het_model_B,
#'                          data = db_model_B,
#'                          G = G_model_B, model = "model_B", estimation = "NLLS",
#'                          hypothesis = "het_r", z = z, start.val = starting)
#'
#' # Store and print results
#' summary(het_model_B_r)
#' het_model_B_r$centrality
#'
#' # Heterogeneous factor (edge -level)
#' z <- as.numeric(as.character(db_model_B$party))
#'
#' # Specify starting values
#' starting <- c(alpha = 0.242486,
#'               beta_gender1 = -0.229895,
#'               beta_party1 = 0.42848,
#'               beta_nchair1 = 3.0959,
#'               phi_within  = 0.396371,
#'               phi_between = 0.414135)
#'
#' # Fit model
#' party_model_B <- net_dep(formula = f_model_B, data = db_model_B,
#'                          G = G_model_B, model = "model_B",
#'                          estimation = "NLLS", hypothesis = "par",
#'                          z = z, start.val = starting)
#'
#' # Store and print results
#' summary(party_model_B)
#' party_model_B$centrality
#' }
#' # WARNING, This toy example is provided only for runtime execution.
#' # Please refer to previous examples for sensible calculations.
#' data("db_alumni_test")
#' data("G_model_A_test")
#' db_model_A <- db_alumni_test
#' G_model_A <- G_model_A_test
#' f_model_A <- formula("les ~ dw")
#' lim_model_A_test <- net_dep(formula = f_model_A, data = db_model_A,
#'                        G = G_model_A, model = "model_A", estimation = "NLLS",
#'                        hypothesis = "lim", start.val = c(alpha = 0.09030594,
#'                                                          beta_dw = 1.21401940,
#'                                                          phi = 1.47140647))
#' summary(lim_model_A_test)
#' @import minpack.lm
#' @importFrom  bbmle parnames mle2
#' @importFrom dplyr %>% group_by summarize
#' @export
net_dep <- function(formula = formula(),
                    data = list(),
                    G = list(),
                    model = c("model_A","model_B"),
                    estimation = c("NLLS","MLE"),
                    hypothesis = c("lim","het","het_l","het_r","par", "par_split_with", "par_split_btw", "par_split_with_btw"),
                    endogeneity = FALSE,
                    correction = NULL,
                    first_step = NULL,
                    z = NULL,
                    formula_first_step = NULL,
                    exclusion_restriction = NULL,
                    start.val = NULL,
                    to_weight = NULL,
                    time_fixed_effect = NULL,
                    ind_fixed_effect = NULL,
                    mle_controls = NULL,
                    kappa = NULL,
                    delta = NULL) {

  if (missing(kappa) & hypothesis != "het_l") {
    kappa = max(rowSums(G))
  } else {
    kappa = max(colSums(G))
  }

  data_list <- prepare_data(formula, data, time_fixed_effect)
  X <- data_list[["data"]]
  y <- data_list[["y"]]
  y_name <- data_list[["y_name"]]
  factors <- data_list[["factors"]]
  tt <- data_list[["time"]]

  if (!is.null(ind_fixed_effect)) {
    n_fe <- length(unique(data[, ind_fixed_effect]))
  } else {
    n_fe <- 0
  }

  env <- environment()
  check_entries(as.list.environment(env))

  e <- toList(X)
  to_rm <- which(names(e) %in% "Ones")

  if (is.null(start.val)) {
    starting.values <- setNames(rep(0.01, length(e)), c("alpha",
                                                        names(e)[ - to_rm]))
    starting.values <- toList(starting.values)
  } else {
    starting.values <- start.val
  }

  n <- nrow(G)
  I <- diag(n)
  e[["I"]] <- I
  e[[y_name]] <- y

  if (endogeneity == TRUE & is.null(first_step)) {
    first_step <- "standard"
  }

  if (endogeneity == TRUE & !is.null(ind_fixed_effect)) {
    first_step <- "fe"
  }

  check_first_step <- !is.null(formula_first_step)

  if (check_first_step) {
    if (!is.null(time_fixed_effect)) {
      formula_first_step_vars <- all.vars(formula_first_step)
      if (time_fixed_effect %in% formula_first_step_vars == F) {
        formula_first_step <- formula(paste0(formula_first_step, " + ", time_fixed_effect))
      }
    }
  }

  if (endogeneity) {

    if (check_first_step) {
      data_first_step_list <- prepare_data(formula_first_step, data, time_fixed_effect)
      X_first_step <- data_first_step_list[["data"]]
      factors_first_step <- data_first_step_list[["factors"]]
      if (is.null(exclusion_restriction)) {
        e_first_step <- create_dyadic_db_no_excl(data = X_first_step[, - 1], G = G,
                                                 option = first_step, tt = tt, is_dummy = factors_first_step)
      } else {
        e_first_step <- create_dyadic_db(data = X_first_step[, - 1], G = G,
                                         exclusion_restriction = exclusion_restriction,
                                         option = first_step, tt = tt, is_dummy = factors_first_step)
      }
    } else {

      if (is.null(exclusion_restriction)) {
        e_first_step <- create_dyadic_db_no_excl(data = X[, - 1], G = G,
                                         option = first_step, tt = tt, is_dummy = factors)
      } else {
        e_first_step <- create_dyadic_db(data = X[, - 1], G = G,
                                         exclusion_restriction = exclusion_restriction,
                                         option = first_step, tt = tt, is_dummy = factors)
      }
    }

    if (first_step == "fe") {
      db_first_step <- e_first_step
      formula_first_step <- formula(paste0(colnames(db_first_step)[1], " ~ ",
                                           paste(colnames(db_first_step)[-1],
                                                 collapse = " + ")))
    } else {
      db_first_step <- e_first_step[, - which(colnames(e_first_step) %in%
                                                c("cidx", "ridx"))]
      formula_first_step <- formula(paste0(colnames(db_first_step)[1], " ~ ",
                                           paste(colnames(db_first_step)[-1],
                                                 collapse = " + ")))
    }

    first_step <- lm(formula = formula_first_step, data = e_first_step)
    res_form <- residuals(first_step)

    if (is.na(sum(res_form))) {
      stop('First step cannot be identified. Try changing specification.')
    }

    if (length(res_form) < nrow(e_first_step)) {
      stop('First step cannot be identified. Try changing specification.')
    }

    if (correction == "heckman") {

      if (!is.null(tt)) {

        if (check_first_step == FALSE) {
          fed <- colnames(X)[which(factors %in% "time_fixed_effect") + 1]
          time <- rep(0, nrow(e_first_step))
        } else {
          fed <- colnames(X_first_step)[which(factors_first_step %in% "time_fixed_effect") + 1]
          time <- rep(0, nrow(e_first_step))
        }

        for (i in 1:length(fed)) {
          time <- time + ifelse (e_first_step[, fed[i]] == 1, i + 1, 0)
        }

        time <- ifelse (time==0, 1, time)
        unobs <- data.frame(cidx = e_first_step$cidx, time = time,
                            unobservables = res_form)
        unobs <- unobs %>% group_by(cidx, time) %>%
          summarize(unobservables = sum(unobservables))
        unobs <- data.frame(unobs)
        colnames(unobs)[1] <- "id"
        to_merge <- data.frame(id = rownames(G), time = tt)
        to_merge <- plyr::join(x = to_merge, y = unobs, by = c("id", "time"), type = "left", match = "all")
        e[["unobservables"]] <- to_merge[, "unobservables"]

      } else {
        unobs <- data.frame(cidx = e_first_step$cidx, unobservables = res_form)
        unobs <- unobs %>% group_by(cidx) %>%
          summarize(unobservables = sum(unobservables))
        unobs <- data.frame(unobs)
        e[["unobservables"]] <- unobs[, "unobservables"]
      }

    } else if (correction=="iv") {
      pred_from <- predict(first_step)
      pred_form_mat <- vec_to_mat(vec = res_form, e_first_step =
                                    e_first_step, fed = NULL, id = rownames(G))
      e[["G"]] <- as.matrix(pred_form_mat)
    }
  }

  if (!is.null(ind_fixed_effect)) {
    individual_fixed_effect <- as.factor(data[, ind_fixed_effect])
    mm <- model.matrix(~individual_fixed_effect - 1, model.frame(~individual_fixed_effect - 1))
    mm <- matrix(mm, nrow(mm), ncol(mm))
    colnames(mm) <- levels(individual_fixed_effect)
    mm <- mm[, -1]
    if (is.null(start.val)) {
      start_val_mm <- setNames(rep(0.01, ncol(mm)), colnames(mm))
      start_val_mm <- toList(start_val_mm)
      starting.values <- c(starting.values, start_val_mm)
    }
    mm_df <- data.frame(mm)
    e <- c(e, as.list(mm_df))
    X <- cbind(X, mm_df)
  }

  if (!is.null(z) & hypothesis != "par" & hypothesis != "par_split_btw" & hypothesis != "par_split_with" & hypothesis != "par_split_with_btw") {
    check_z <- Reduce("+", lapply(X, identical, z))
    if (check_z == 0) {
      X <- cbind(X, z)
      e[["z"]] <- z
      if (is.null(start.val)) {
        starting.values[["z"]] <- 0.01
      }
    }
  }

  set_second_step <- nls_prepare_data(dependent_variable = y_name, X = X,
                                      G = G, z = z, e = e,
                                      start.val = start.val,
                                      starting.values = starting.values,
                                      model = model, hypothesis = hypothesis,
                                      endogeneity = endogeneity,
                                      correction = correction, tt = tt)
  nls_formula_fit <- set_second_step[[1]]
  starting.values <- set_second_step[[2]]

  if (estimation == "NLLS" & is.null(delta)) {

    e <- set_second_step[[3]]
    set.control <- nls.lm.control(ftol = sqrt(.Machine$double.eps) / 1000,
                                  ptol = sqrt(.Machine$double.eps) / 1000,
                                  gtol = 0, diag = list(), epsfcn = 0,
                                  factor = 100, maxfev = integer(),
                                  maxiter = 300, nprint = 0)

    if (!is.null(to_weight)) {
      e[["to_weight"]] <- to_weight
      second_step <- nlsLM(nls_formula_fit, start = starting.values, data = e,
                           trace = T, weights = to_weight,
                           control = set.control)
      second_step$data <- e
      second_step$m$mformula <- NULL
      second_step$delta <- FALSE
    } else {
      second_step <- nlsLM(nls_formula_fit, start = starting.values, data = e,
                           trace= T, control = set.control)
      second_step$data <- e
      second_step$m$mformula <- NULL
      second_step$delta <- FALSE
    }

  } else if (estimation == "NLLS" & !is.null(delta)) {

    e <- set_second_step[[3]]
    set.control <- nls.lm.control(ftol = sqrt(.Machine$double.eps) / 1000,
                                  ptol = sqrt(.Machine$double.eps) / 1000,
                                  gtol = 0, diag = list(), epsfcn = 0,
                                  factor = 100, maxfev = integer(),
                                  maxiter = 300, nprint = 0)

    nls_const <- nls_constrained(nls_formula_fit = nls_formula_fit, G = G, e = e,
                    hypothesis = hypothesis, model = model, delta = delta, to_weight = to_weight)

    nls_formula_fit_bou <- nls_const$nls_formula_fit
    e_new <- nls_const$e

    if (!is.null(to_weight)) {
      e_new[["to_weight"]] <- nls_const$to_weight
      second_step <- nlsLM(nls_formula_fit_bou, start = starting.values, data = e_new,
                           trace = T, weights = to_weight,
                           control = set.control)
      second_step$data <- e
      second_step$m$mformula <- nls_formula_fit
      second_step$delta <- TRUE
    } else {
      second_step <- nlsLM(nls_formula_fit_bou, start = starting.values, data = e_new,
                           trace= T, control = set.control)
      second_step$data <- e
      second_step$m$mformula <- nls_formula_fit
      second_step$delta <- TRUE
    }

  } else if (estimation == "MLE") {
    set_second_step <- mle_prepare_data(X = X, G = G, z = z, e = e,
                                        start.val = start.val,
                                        starting.values = starting.values,
                                        model = model, hypothesis = hypothesis,
                                        endogeneity = endogeneity,
                                        correction = correction, tt = tt,
                                        mle_controls = mle_controls,
                                        kappa = kappa)
    starting.values <- set_second_step[[1]]
    boundU <- set_second_step[[2]]
    boundL <- set_second_step[[3]]
    e <- set_second_step[[4]]
    if (!is.null(ind_fixed_effect)) {
      X <- X[, -1]
    }

    if (endogeneity) {
      X[, "unobservables"] <- e[["unobservables"]]
    }

    if (model == "model_A") {
      if (hypothesis == "lim") {
        second_step <- mle_compet_lim(y = y, X = X, G = G, z = z,
                                      starting.values = starting.values,
                                      boundL = boundL, boundU = boundU, ind_fixed_effect)
      } else if (hypothesis == "het") {
        second_step <- mle_compet_het(y = y, X = X, G = G, z = z,
                                      starting.values = starting.values,
                                      boundL = boundL, boundU = boundU)
      }
    }

    if (model == "model_B") {
      if (hypothesis == "lim") {
        second_step <- mle_compet_lim(y = y, X = X, G = G, z = z,
                                      starting.values = starting.values,
                                      boundL = boundL, boundU = boundU, ind_fixed_effect)
      } else if (hypothesis == "het_l") {
        second_step <- mle_het(y = y, X = X, G = G, z = z, side = "het_l",
                               starting.values = starting.values,
                               boundL = boundL, boundU = boundU)
      } else if (hypothesis == "het_r") {
        second_step <- mle_het(y = y, X = X, G = G, z = z, side = "het_r",
                               starting.values = starting.values,
                               boundL = boundL, boundU = boundU)
      } else if (hypothesis == "par") {
        second_step <- mle_party(y = y, X = X, G_within = e[["G_within"]],
                                 G_between = e[["G_between"]],
                                 starting.values = starting.values,
                                 boundL = boundL, boundU = boundU)
      } else if (hypothesis == "par_split_with") {
        second_step <- mle_party_split_with(y = y, X = X,
                                            G_within_0 = e[["G_within_0"]],
                                            G_within_1 = e[["G_within_1"]],
                                            G_between = e[["G_between"]],
                                            starting.values = starting.values,
                                            boundL = boundL, boundU = boundU)
      } else if (hypothesis == "par_split_btw") {
        second_step <- mle_party_split_btw(y = y, X = X,
                                           G_within = e[["G_within"]],
                                           G_between_01 = e[["G_between_01"]],
                                           G_between_10 = e[["G_between_10"]],
                                           starting.values = starting.values,
                                           boundL = boundL, boundU = boundU)
      } else if (hypothesis == "par_split_with_btw") {
        second_step <- mle_party_split_with_btw(y = y, X = X,
                                                G_within_0 = e[["G_within_0"]],
                                                G_within_1 = e[["G_within_1"]],
                                                G_between_01 = e[["G_between_01"]],
                                                G_between_10 = e[["G_between_10"]],
                                                starting.values = starting.values,
                                                boundL = boundL, boundU = boundU)
      }
    }

    second_step@formula <- paste(nls_formula_fit[2], nls_formula_fit[3], sep=' ~ ')

  }

   centrality <- parameter_dependent_centrality(second_step = second_step,
                                                hypothesis = hypothesis, I = I,
                                                G = G, e = e)

  if (endogeneity == TRUE) {
    res <- list(second_step = second_step, centrality = centrality,
                first_step = first_step, hypothesis = hypothesis,
                starting.values = starting.values)
  } else {
    res <- list(second_step = second_step, centrality = centrality,
                first_step = NULL, hypothesis = hypothesis,
                starting.values = starting.values)
  }

  class(res) <- "econet"
  attributes(res)$hypothesis <- hypothesis
  attributes(res)$model <- model
  return(res)

}
