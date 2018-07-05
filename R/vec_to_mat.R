#' Function to transform a vector in a matrix
#'
#' @param vec a numeric vector.
#' @param e_first_step an object of class \code{data.frame}. It contains the dyadic version of \code{data} provided in \code{net_dep}.
#' @param fed a vector of positive integers. It is used to index fixed effect variables in \code{e_first_step}.
#' @param id a vector of positive integers used as ids for the agents.
#' @return an object of class \code{matrix}.
#' @keywords internal
#' @import Matrix
#' @importFrom igraph E
#' @importFrom Matrix bdiag
#' @noRd
vec_to_mat <- function(vec, e_first_step, fed, id) {

  if (!is.null(fed)) {
    n <- e_first_step[, fed]
    baseline <- c(apply(n, 1, sum))
    baseline <- ifelse(baseline == 0, 1, 0)
    n <- cbind(baseline, n)
    mm_list <- lapply(n, function(y) {
      sel <- which(y == 1)
      tmp_id <- id[sel]
      df <- data.frame(u = vec[sel], i = e_first_step[sel, "cidx"],
                       j = e_first_step[sel, "ridx"], stringsAsFactors = F)
      df[,"i"] <- as.character(df[, "i"])
      df[,"j"] <- as.character(df[, "j"])
      g <- igraph::graph_from_edgelist(el = as.matrix(df[, c("i", "j")]),
                                       directed = T)
      E(g)$weight <- df[, "u"]
      mm <- igraph::as_adjacency_matrix(g, type = "both", attr = "weight")
      ord <- order(match(rownames(mm), rownames(tmp_id)))
      mm <- mm[ord, ord]
      mm
    })
    mm <- bdiag(mm_list)
  } else {
    df <- data.frame(u = vec, i = e_first_step[,"cidx"],
                     j = e_first_step[, "ridx"], stringsAsFactors = F)
    df[,"i"] <- as.character(df[, "i"])
    df[,"j"] <- as.character(df[, "j"])
    g <- igraph::graph_from_edgelist(el = as.matrix(df[,c("i","j")]),
                                     directed = T)
    E(g)$weight <- df[,"u"]
    mm <- igraph::as_adjacency_matrix(g, type = "both", attr = "weight")
  }

  return(as.matrix(mm))

}
