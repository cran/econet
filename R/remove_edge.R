#' Remove one edge in the adjacency matrix
#'
#' @param x string. id of agent \emph{i}
#' @param y string. id of agent \emph{j}
#' @param mat an object of class \code{matrix}: adjacency matrix.
#' @return return an object of class \code{igraph} where the edge connecting \emph{i} and \emph{j} in the network is removed.
#' @keywords internal
#' @importFrom igraph graph.adjacency
#' @noRd
remove_edge <- function(x, y, mat){
  i_x <- which(colnames(mat) %in% x)
  i_y <- which(colnames(mat) %in% y)
  tmp_mat <- mat
  tmp_mat[i_x, i_y] <- tmp_mat[i_y, i_x] <- 0
  graph.adjacency(ifelse(tmp_mat > 0, 1, 0), mode = "undirected")
}
