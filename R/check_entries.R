#' Check if  the entries in \code{net_dep} are correctly provided by the user.
#'
#' @param env list of objects passed through the arguments of \code{net_dep}
#' @return It returns an error message if  entries are incorrectly used.
#' @keywords internal
#' @noRd
check_entries <- function(env) {

  for (i in 2:length(env)) {
    assign(names(env)[i], env[[i]], environment())
  }

  if (((model %in% c("hypothesis_A", "hypothesis_B")) |
      (estimation %in% c("NLLS", "MLE")) |
      (hypothesis %in% c("lim", "het", "het_l", "het_r", "par")) |
      (as.character(endogeneity) %in% c("TRUE", "FALSE"))) == FALSE) {
    stop('Incorrect field defined')
  }

  if (!is.null(correction)) {
    if (correction %in% c("heckman", "iv") == FALSE) {
      stop('Incorrect field defined')
    }
  }

  if (!is.null(first_step)) {
    if (first_step %in% c("standard", "fe", "shortest", "coauthors", "degree") == FALSE) {
      stop('Incorrect field defined')
    }
  }

  if (is.null(data)) {
    stop('Data is missing')
  }

  if (is.null(G)) {
    stop('G is missing')
  }

  if (nrow(G) != ncol(G)) {
    stop('G must be squared')
  }

  if (!is.null(z) & !is.vector(z)) {
    stop('z must be a vector')
  }

  if (!is.null(z) & !is.numeric(z)) {
    stop('z must be numeric')
  }

  if (!is.null(to_weight) & !is.vector(to_weight)) {
    stop('to_weight must be a vector')
  }

  if (!is.null(to_weight) & !is.numeric(to_weight)) {
    stop('to_weight must be numeric')
  }

  if (!is.null(tt) & !is.vector(tt)) {
    stop('tt must be a vector')
  }

  if (!is.null(tt) & !is.vector(tt)) {
    stop('tt must be numeric')
  }

  if (sum(hypothesis%in%c("het", "het_l", "het_r", "par")) > 0 & is.null(z)) {
    stop('z is missing')
  }

  if (!is.null(first_step) & is.null(exclusion_restriction)) {
    stop('exclusion_restriction is missing')
  }

  if (is.na(sum(y + X + G)) |
     (!is.null(z) & is.na(sum(z))) |
     (!is.null(tt) & is.na(sum(tt))) |
     (!is.null(to_weight) & is.na(sum(to_weight))) |
     (!is.null(data_first_step) & is.na(sum(data_first_step))) |
     (!is.null(start.val) & is.na(sum(Reduce("c",start.val)))) |
     (!is.null(exclusion_restriction) & is.na(sum(exclusion_restriction)))) {
    stop('net_dep does not know how to deal with missing values')
  }

  if (!is.null(first_step) ) {
    if  (is.null(rownames(G))) {
      stop('row and column names in G cannot be null. Hint: they should indicate nodes id')
    }

  }

  if (nrow(data) != nrow(G)) {
    stop('Data are not conformable')
  }

  if (hypothesis != "lim" & !is.null(z)) {
    if (hypothesis != "lim" & length(z) != nrow(G)) {
      stop('Data are not conformable')
    }
  }

  if (!is.null(exclusion_restriction)) {
    if (nrow(exclusion_restriction) != nrow(G)) {
      stop('Data are not conformable')
    }
  }

  if (!is.null(tt)) {
    if (length(tt) != nrow(G)) {
      stop('Data are not conformable')
    }
  }

  if (!is.null(to_weight)) {
    if (length(to_weight) != nrow(G)) {
      stop('Data are not conformable')
    }
  }

  if ((model == "hypothesis_A" & sum(hypothesis%in%c("het_l", "het_r")) > 0) |
     (model == "hypothesis_B" & hypothesis=="het")) {
    stop('model is incompatible with hypothesis')
  }

  if (endogeneity == FALSE) {
    parameters <- c("alpha", "phi", "gamma", "theta_0", "theta_1", "eta_0",
                    "eta_1", "phi_within", "phi_between")
    if (
      (!is.null(start.val) & hypothesis == "lim" &
       sum(c("alpha", "phi") %in% names(start.val)) != 2 &
       sum(parameters %in% names(start.val)) != 2 ) |
      (!is.null(start.val) & hypothesis == "het" &
       sum(c("alpha", "phi", "gamma") %in% names(start.val)) != 3 &
       sum(parameters %in% names(start.val)) != 3 ) |
      (!is.null(start.val) & hypothesis == "het_l" &
       sum(c("alpha", "theta_0", "theta_1")%in%names(start.val)) != 3 &
       sum(parameters %in% names(start.val)) != 3 ) |
      (!is.null(start.val) & hypothesis == "het_r" &
       sum(c("alpha", "eta_0","eta_1")
        %in%names(start.val)) != 3 &
       sum(parameters %in% names(start.val)) != 3 ) |
      (!is.null(start.val) & hypothesis == "par" &
       sum(c("alpha", "phi_within", "phi_between") %in% names(start.val)) != 3 &
       sum(parameters %in% names(start.val)) != 3)) {
      stop('Starting values are not correctly specified')
    }

    if ((hypothesis=="lim" & estimation == "NLLS" & !is.null(start.val) &
         length(start.val) > (ncol(X)+2)) |
       (hypothesis=="lim" & estimation == "MLE" & !is.null(start.val) &
        length(start.val) > (ncol(X)+3))) {
      stop('Starting values are not correctly specified')
    }

    if ((hypothesis != "lim" & estimation == "NLLS" & !is.null(start.val) &
        length(start.val) > (ncol(X) + 3)) |
       (hypothesis != "lim" & estimation == "MLE" & !is.null(start.val) &
        length(start.val) > (ncol(X) + 4))) {
      stop('Starting values are not correctly specified')
    }

  } else {
    parameters <- c("alpha", "phi", "gamma", "theta_0", "theta_1", "eta_0",
                    "eta_1", "phi_within", "phi_between", "unobservables")
    if (
      (!is.null(start.val) & hypothesis == "lim" &
       sum(c("alpha", "phi") %in% names(start.val)) != 2 &
       sum(parameters %in% names(start.val)) != 3 ) |
      (!is.null(start.val) & hypothesis == "het" &
       sum(c("alpha", "phi", "gamma") %in% names(start.val)) != 3 &
       sum(parameters %in% names(start.val)) != 4 ) |
      (!is.null(start.val) & hypothesis == "het_l" &
       sum(c("alpha", "theta_0", "theta_1") %in% names(start.val)) != 3 &
       sum(parameters %in% names(start.val)) != 4 ) |
      (!is.null(start.val) & hypothesis == "het_r" &
       sum(c("alpha", "eta_0", "eta_1") %in% names(start.val)) != 3 &
       sum(parameters %in% names(start.val)) != 4 ) |
      (!is.null(start.val) & hypothesis == "par" &
       sum(c("alpha", "phi_within", "phi_between") %in% names(start.val)) != 3 &
       sum(parameters %in% names(start.val)) != 4)) {
      stop('Starting values are not correctly specified')
    }

    if ((hypothesis=="lim" & estimation == "NLLS" & !is.null(start.val) &
         length(start.val) > (ncol(X) + 3)) |
       (hypothesis=="lim" & estimation == "MLE" & !is.null(start.val) &
        length(start.val) > (ncol(X) + 4))) {
      stop('Starting values are not correctly specified')
    }

    if ((hypothesis != "lim" & estimation == "NLLS" & !is.null(start.val) &
         length(start.val) > (ncol(X) + 4)) |
       (hypothesis != "lim" & estimation == "MLE" & !is.null(start.val) &
        length(start.val) > (ncol(X) + 5))) {
      stop('Starting values are not correctly specified')
    }

  }

  if (endogeneity == TRUE & is.null(correction)) {
    stop('if  endogeneity == TRUE, correction must be specified')
  }

  if (!is.null(mle_controls) & is.null(names(mle_controls))) {
    stop('mle_controls is not correctly specified (names must be provided)')
  }

}

