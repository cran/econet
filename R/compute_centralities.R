#' Compute traditional measures of network centrality.
#'
#' @param mat an object of class \code{matrix}: adjacency matrix.
#' @param directed logical. If it is \code{TRUE} the network is considered directed.
#' @param weighted logical. If it is \code{TRUE} the network is considered weighted.
#' @param normalization string. Default is FALSE. "bygraph" and "bycomponent" divide degree and closeness by \eqn{n-1} and betweenness by \eqn{(n-1)*(n-2)}, where n is equal to the number of nodes in the whole network, or in the component. "bymaxgraph" and "bymaxcomponent" divide the vector of centrality measures by the maximum value of centrality in the whole network or in the component.
#' @return an object of class \code{data.frame} containing the vectors of centrality measures.
#' @keywords internal
#' @importFrom intergraph asNetwork
#' @importFrom sna evcent
#' @importFrom igraph E V
#' @noRd
compute_centralities <- function(mat,
                                 directed = c("TRUE","FALSE"),
                                 weighted = FALSE,
                                 normalization = FALSE) {

  direction <- ifelse(directed == TRUE, "directed", "undirected")
  weighted <- ifelse(is.null(weighted), FALSE, weighted)
  normalization <- ifelse(is.null(normalization), FALSE, normalization)

  if (weighted == FALSE) {
    g <- igraph::graph_from_adjacency_matrix(adjmatrix = ifelse(mat > 0, 1, 0),
                                             mode = direction)
  } else {
    mat <- 1/mat
    mat[is.infinite(mat)] <- 0
    g <- igraph::graph_from_adjacency_matrix(adjmatrix = mat, mode = direction,
                                             weighted = TRUE)
  }

  if (directed == FALSE) {
    indegree <- outdegree <- igraph::degree(graph = g)
    betweenness <- igraph::betweenness(graph = g, directed = FALSE)
    eigenvector <- sna::evcent(asNetwork(g))

    isolate <- c(which(indegree == 0 | outdegree == 0))
    if (length(isolate)>0) {
      outcloseness <- rep(0, igraph::vcount(g))
      incloseness <- rep(0, igraph::vcount(g))
      g <- g - V(g)[V(g)$name%in%names(isolate)]
      outcloseness[-as.numeric(isolate)] <- suppressWarnings(igraph::closeness(graph = g,
                                                              mode = "all"))
      incloseness[-as.numeric(isolate)] <- suppressWarnings(igraph::closeness(graph = g,
                                                             mode = "all"))
    } else {
      outcloseness <- suppressWarnings(igraph::closeness(graph = g, mode = "all"))
      incloseness <- suppressWarnings(igraph::closeness(graph = g, mode = "all"))
    }
  } else {

    indegree <- igraph::degree(graph = g, mode = "in")
    outdegree <- igraph::degree(graph = g, mode = "out")
    betweenness <- igraph::betweenness(graph = g)
    eigenvector <- sna::evcent(asNetwork(g))

    isolate <- c(which(indegree == 0 | outdegree == 0))
    if (length(isolate)>0) {
      outcloseness <- rep(0, igraph::vcount(g))
      incloseness <- rep(0, igraph::vcount(g))
      g <- g - V(g)[V(g)$name%in%names(isolate)]
      outcloseness[-as.numeric(isolate)] <- suppressWarnings(igraph::closeness(graph = g,
                                                              mode = "out"))
      incloseness[-as.numeric(isolate)] <- suppressWarnings(igraph::closeness(graph = g,
                                                             mode = "in"))
    } else {
      outcloseness <- suppressWarnings(igraph::closeness(graph = g, mode = "out"))
      incloseness <- suppressWarnings(igraph::closeness(graph = g, mode = "in"))
    }
  }

  res <- data.frame(gov_id = rownames(mat), indegree = indegree,
                    outdegree = outdegree, betweenness = betweenness,
                    incloseness = incloseness, outcloseness = outcloseness,
                    eigenvector = eigenvector, stringsAsFactors = F)

    if (normalization == "bycomponent" | normalization == "bymaxcomponent") {
      c_g <- igraph::components(g)$membership
      name_c_g <- unique(as.numeric(c_g))
      for(i in 1:length(name_c_g)) {
        id_sel <- names(c_g)[c_g %in% name_c_g[i]]
        res_sel <- which(res[,"gov_id"] %in% id_sel)
        n <- length(id_sel)
        if (n>2 & normalization == "bycomponent") {
          res[res_sel, "indegree"] <- res[res_sel, "indegree"] / (n - 1)
          res[res_sel, "outdegree"] <- res[res_sel, "outdegree"] / (n - 1)
          res[res_sel, "outcloseness"] <- res[res_sel, "outcloseness"] / (n - 1)
          res[res_sel, "incloseness"] <- res[res_sel, "incloseness"] / (n - 1)
          if (directed == FALSE) {
            res[res_sel, "betweenness"] <- res[res_sel, "betweenness"] /
              (((n - 1)*(n - 2)) / 2)
          } else {
            res[res_sel,"betweenness"] <- res[res_sel,"betweenness"] /
              (((n - 1) * (n - 2)))
          }
        } else if (n > 1 & normalization == "bymaxcomponent") {
          res[res_sel, "indegree"] <- res[res_sel, "indegree"] /
            max(res[res_sel, "indegree"])
          res[res_sel, "outdegree"] <- res[res_sel, "outdegree"] /
            max(res[res_sel, "outdegree"])
          res[res_sel, "betweenness"] <- res[res_sel, "betweenness"] /
            max(res[res_sel, "betweenness"])
          res[res_sel, "outcloseness"] <- res[res_sel, "outcloseness"] /
            max(res[res_sel, "outcloseness"])
          res[res_sel, "incloseness"] <- res[res_sel," incloseness"] /
            max(res[res_sel, "incloseness"])
        }
      }
    } else if (normalization == "bygraph") {
      n <- nrow(res)
      res[, "indegree"] <- res[, "indegree"] / (n - 1)
      res[, "outdegree"] <- res[, "outdegree"]/(n - 1)
      res[, "outcloseness"] <- res[, "outcloseness"] / (n - 1)
      res[, "incloseness"] <- res[, "incloseness"] / (n - 1)
      if (directed == FALSE) {
        res[, "betweenness"] <- res[, "betweenness"] /
          (((n - 1) * (n - 2)) / 2)
      } else {
        res[, "betweenness"] <- res[, "betweenness"] / (((n - 1) * (n - 2)))
      }
    } else if (normalization == "bymaxgraph") {
      res[, "indegree"] <- res[, "indegree"] / max(res[, "indegree"])
      res[, "outdegree"] <- res[, "outdegree"] / max(res[, "outdegree"])
      res[, "betweenness"] <- res[, "betweenness"] / max(res[, "betweenness"])
      res[, "outcloseness"] <- res[, "outcloseness"] / max(res[, "outcloseness"])
      res[, "incloseness"] <- res[, "incloseness"] / max(res[, "incloseness"])
    }

  res[which(is.na(res[, "indegree"])), "indegree"] <- 0
  res[which(is.na(res[, "outdegree"])), "outdegree"] <- 0
  res[which(is.na(res[, "betweenness"])), "betweenness"] <- 0
  res[which(is.na(res[, "outcloseness"])), "outcloseness"] <- 0
  res[which(is.na(res[, "incloseness"])), "incloseness"] <- 0
  res[which(is.na(res[, "eigenvector"])), "eigenvector"] <- 0
  rownames(res) <- 1:nrow(res)
  return(res)
}
