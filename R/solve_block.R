#' Invert a block diagonal matrix block by block
#'
#' @param mat an object of class \code{matrix}: adjacency matrix.
#' @return an object of class \code{matrix}: inverted adjacency matrix.
#' @keywords internal
#' @importFrom MASS ginv
#' @noRd
solve_block <- function(mat){

    x <- mat
    diag(x) <- 1
    g <- igraph::graph_from_adjacency_matrix(adjmatrix = as.matrix(x),
                  mode = "directed", weighted = TRUE, diag = TRUE)
    groups <- Map(sort, igraph::neighborhood(g, nrow(mat)))
    groups <- unique(Map(as.numeric, groups))
    sub.Mat <- Map(`[`, list(mat), groups, groups, drop = FALSE)
    inv.sub.Mat <- Map(ginv,sub.Mat)
    mat.inv <- matrix(0, nrow(mat), ncol(mat))
    for (i in 1:length(sub.Mat)) {
      mat.inv[groups[[i]], groups[[i]]] <- inv.sub.Mat[[i]]
    }

  return(mat.inv)

}
