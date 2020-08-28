#' Prepare \code{data} to be elaborated by \code{net_dep}.
#'
#' @param formula an object of class \code{formula}: a symbolic description of the model to be fitted. The constant (i.e. intercept) and the autogressive parameter needs not to be specified.
#' @param time_fixed_effect an optional string. It indicates the name of the time index used in formula. It is used for models with longitudinal data.
#' @keywords internal
#' @importFrom stats model.frame model.response model.matrix
#' @noRd
prepare_data <- function(formula, data, time_fixed_effect){

  m <- model.frame(formula, data = data)

  name_controls <- all.vars(formula)[-1]

  y <- model.response(data = m, type = "numeric")
  y_name <- colnames(m)[1]

  X <- model.matrix(m, data = data)
  colnames(X) <- c("Ones", colnames(X)[ - 1])

  contrast_names <- names(attr(X, "contrasts"))

  factors <- rep(name_controls, table(attr(X, "assign")[ - 1]))
  factors[factors %in% time_fixed_effect] <- "time_fixed_effect"
  factors[factors %in% contrast_names] <- "factor"
  factors[-which(factors %in% c("time_fixed_effect","factor"))] <- "numeric"

  X <- data.frame(X)

  if (!is.null(time_fixed_effect)) {
    tt <- data[, time_fixed_effect]
    tt <- as.factor(as.numeric(tt))
    tt <- as.numeric(as.character(tt))
  } else {
    tt <- NULL
  }

  prepared_data <- list()
  prepared_data[["data"]] <- X
  prepared_data[["y"]] <- y
  prepared_data[["y_name"]] <- y_name
  prepared_data[["factors"]] <- factors
  prepared_data[["time"]] <- tt

  return(prepared_data)
}


