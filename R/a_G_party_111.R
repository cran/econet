#' Column-normalized adjacency matrix representing the alumni network of the 111th U.S. Congress weighted by party affiliation. Weights take values of 1 if the linked congressmen have the same school and party, 0.5 if they have the same school and a different party, 0 if they have a different school and the same party, and -0.5 if they have a different school and a different party.
#'
#' @docType data
#' @usage data("a_G_party_111")
#' @format an object of class \code{Matrix} with 426 rows and 426 columns
#' @keywords datasets
#' @references
#' Battaglini M., E. Patacchini (2018), "Influencing Connected Legislators", Journal of Political Economy, forthcoming.\cr
"a_G_party_111"

