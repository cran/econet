#' Create a dyadic version of \code{data} provided by the user in \code{net_dep}.
#'
#' @param data an object of class \code{data.frame} containing the list of control variables specified by the user.
#' @param G an object of class \code{Matrix} representing the social network.
#' @param exclusion_restriction an object of class \code{Matrix} representing the exogenous matrix used to instrument the endogenous social network, if \code{endogeneity = TRUE}.
#' @param option string. It contains the model specification indicated by the user in the argument \code{hypothesis} through \code{net_dep}.
#' @param tt a vector of positive integers used as time index. It is used for models with longitudinal data.
#' @param is_dummy a vector of integers used to index which variables in \code{data} are factors.
#' @return an object of class \code{data.frame} containing a dyadic version of \code{data}.
#' @keywords internal
#' @noRd
create_dyadic_db <- function(data, G, exclusion_restriction,
                             option = c("standard", "fe", "shortest",
                                        "coauthors", "degree"), tt, is_dummy) {

  x_names <- colnames(data)

  if (!is.null(tt)) {
    k <- length(unique(tt))
    z <- ncol(data)
    list_net_form_data <- list()

    for (i in 1:k) {
      tmp_data <- data[tt == i, ]
      tmp_el_g <- mat_to_el(G[tt == i, tt == i])
      tmp_el_w <- mat_to_el(exclusion_restriction[tt == i, tt == i])
      tmp_df <- matrix("", nrow = length(tmp_el_g), ncol = z)
      tmp_df <- data.frame(tmp_df)
      colnames(tmp_df) <- x_names

      for (n in 1:z) {
        x <- tmp_data[, colnames(tmp_data) %in% x_names[n]]
        if (is_dummy[n] == "numeric") {
          x <- as.numeric(as.character(x))
          tmp_mat_x <- abs(outer(x, x, "-"))
          tmp_mat_x <- mat_to_el(tmp_mat_x)
          tmp_df[,n] <- tmp_mat_x
        } else if (is_dummy[n] == "factor") {
          tmp_mat_x <- outer(x, x, "ifdummy")
          tmp_mat_x <- mat_to_el(tmp_mat_x)
          tmp_df[,n] <- tmp_mat_x
        } else if (is_dummy[n] == "time_fixed_effect") {
          tmp_mat_x <- outer(x, x, "ifdummy1")
          tmp_mat_x <- mat_to_el(tmp_mat_x)
          tmp_df[,n] <- tmp_mat_x
        }
      }

      res <- data.frame(y = tmp_el_g, exclusion_restriction = tmp_el_w, tmp_df)
      list_net_form_data[[i]] <- res
    }

    data_network_formation <- Reduce("rbind", list_net_form_data)
    fixed <- Reduce("rbind", lapply(lapply(1:k, function(y)
      rownames(G)[which(tt == y)]), create_dyadic_fe))
    net_form_data <- data.frame(data_network_formation, cidx = fixed[, 1],
                                ridx = fixed[, 2])

  } else {

    z <- ncol(data)
    tmp_data <- data
    tmp_el_g <- mat_to_el(G)
    tmp_el_w <- mat_to_el(exclusion_restriction)
    tmp_df <- matrix("", nrow = length(tmp_el_g), ncol = z)
    tmp_df <- data.frame(tmp_df)
    colnames(tmp_df) <- x_names

    for (n in 1:z) {
      x <- tmp_data[,colnames(tmp_data) %in% x_names[n]]
      if (is_dummy[n] == "numeric") {
        x <- as.numeric(as.character(x))
        tmp_mat_x <- abs(outer(x, x, "-"))
        tmp_mat_x <- mat_to_el(tmp_mat_x)
        tmp_df[,n] <- tmp_mat_x
      } else {
        tmp_mat_x <- outer(x, x, "ifdummy")
        tmp_mat_x <- mat_to_el(tmp_mat_x)
        tmp_df[,n] <- tmp_mat_x
      }
    }

    data_network_formation <- data.frame(y = tmp_el_g, exclusion_restriction =
                                           tmp_el_w, tmp_df)
    fixed <- create_dyadic_fe(rownames(G))
    net_form_data <- data.frame(data_network_formation, cidx = fixed[, 1],
                                ridx = fixed[, 2])

  }

  if (option == "shortest" | option == "coauthors") {

    if (!is.null(tt)) {
      shortest <- list()
      coauthors <- list()

      for (i in 1:k) {

        tmp_mat <- G[tt == i, tt == i]
        tmp_g <- igraph::graph.adjacency(ifelse(tmp_mat > 0, 1, 0),
                                         mode = "undirected")
        el <- data.frame(sender = c(apply(as.matrix(colnames(tmp_mat)), 1,
                                          function(x) rep(x, nrow(tmp_mat)))),
                         receiver = rep(colnames(tmp_mat), nrow(tmp_mat)),
                         stringsAsFactors = FALSE)[lower.tri(tmp_mat), ]
        el <- as.matrix(el)

        if (option == "shortest") {
          shortest[[i]] <- apply(el, 1, function(x)
            igraph::shortest.paths(remove_edge(x[1], x[2], tmp_mat), v = x[1],
                                   to = x[2]))
        } else if (option == "coauthors") {
          coauthors[[i]] <- apply(el, 1, function(x)
            igraph::cocitation(tmp_g, V(tmp_g)[x[1]])[V(tmp_g)[x[2]]])
        }

      }

      if (option == "shortest") {
        shortest <- as.numeric(as.character(Reduce("c", shortest)))
        shortest[is.na(shortest)] <- 0
      } else if (option == "coauthors") {
        coauthors <- as.numeric(as.character(Reduce("c", coauthors)))
        coauthors[is.na(coauthors)] <- 0
      }

    } else {

      tmp_mat <- G
      tmp_g <- igraph::graph.adjacency(ifelse(tmp_mat>0, 1, 0),
                                       mode = "undirected")
      el <- data.frame(sender = c(apply(as.matrix(colnames(tmp_mat)), 1,
                                        function(x) rep(x, nrow(tmp_mat)))),
                       receiver = rep(colnames(tmp_mat), nrow(tmp_mat)),
                       stringsAsFactors = FALSE)[lower.tri(tmp_mat), ]
      el <- as.matrix(el)

      if (option == "shortest") {
        shortest <- apply(el,1, function(x)
          igraph::shortest.paths(remove_edge(x[1], x[2], tmp_mat), v = x[1],
                                 to = x[2]))
        shortest[is.na(shortest)] <- 0
      } else if (option == "coauthors") {
        coauthors <- apply(el, 1, function(x)
          igraph::cocitation(tmp_g, V(tmp_g)[x[1]])[V(tmp_g)[x[2]]])
        coauthors[is.na(coauthors)] <- 0
      }
    }

    if (option == "shortest") {
      net_form_data <- cbind(net_form_data, shortest)
    } else if (option == "coauthors") {
      net_form_data <- cbind(net_form_data, coauthors)
    }
  }

  if (option == "degree") {

    if (!is.null(tt)) {

    deg <- lapply(1:k, function(y) {
        tmp_g <- G[tt == y, tt == y]
        tmp_g <- ifelse(tmp_g > 0, 1, 0)
        tmp_deg <- apply(tmp_g, 1, sum)
        tmp_mat_deg <- abs(outer(tmp_deg, tmp_deg, "-"))
        tmp_mat_deg <- mat_to_el(tmp_mat_deg)
        tmp_mat_deg
      })
    degree <- Reduce("c", deg)
    net_form_data <- cbind(net_form_data, degree)

    } else {
      tmp_g <- ifelse(tmp_g > 0, 1, 0)
      tmp_deg <- apply(tmp_g, 1, sum)
      tmp_mat_deg <- abs(outer(tmp_deg, tmp_deg, "-"))
      degree <- mat_to_el(tmp_mat_deg)
      net_form_data <- cbind(net_form_data, degree)
    }
  }

  return(net_form_data)

}
