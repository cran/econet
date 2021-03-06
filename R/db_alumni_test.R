#' Dataset from Battaglini, Patacchini (2018)
#'
#' @docType data
#' @usage data("db_alumni_test")
#' @format Data set for runtime executions. An object of class \code{data.frame} with 42 rows and 48 columns
#' \describe{
#' \item{id}{unique congressman id.}
#' \item{PAC}{PAC contributions to a member of Congress (in US dollars).}
#' \item{isolate}{dummy variable taking value of 1 if the congressman did not graduate from the same institution within eight years with any other congressman.}
#' \item{S1-S65}{list of school dummies. Each one refers to a different school. The dummy takes value one if the congressman attended the given school.}
#' \item{S_other}{dummy variable taking value of 1 if the congressman attended any other school not previously listed.}
#' \item{party}{dummy variable taking value 1 if the congressman is a Democrat, and 0 if the congressman is a Republican.}
#' \item{gender}{dummy variable taking value of 1 if the congressman is female.}
#' \item{nchair}{dummy variable taking value of 1 if the congressman is a chair of at least one committee.}
#' \item{sen}{number of consecutive years in the House of Representatives.}
#' \item{margin_b}{dummy variable taking value of 1 if the election margin of victory is less than 5\%.}
#' \item{dw}{distance to the center in terms of ideology measured using the absolute value of the first dimension of the dw-nominate score created by McCarty et al. (1997).}
#' \item{majority}{dummy variable taking value of 1 if the member of Congress is a member of the party holding the majority of the seats in the House of Representatives.}
#' \item{relcom}{dummy variable taking value of 1 if the congressman is member of a powerful committee (Appropriations, Budget, Rules and Ways and Means).}
#' \item{time}{categorical variable. It takes 1 if the record refers to the 109th congress, 2 if the record refers to the 110th congress, 3 if the record refers to the 111th congress, 4 if the record refers to the 112th congress, 5 if the record refers to the 113th congress.}
#' \item{weights}{for each congressman \emph{i}, the value is equal to the inverse of the variance of PAC contributions received by all congressmen in \emph{i}'s State of election.}
#' }
#' @keywords datasets
"db_alumni_test"
