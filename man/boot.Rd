% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boot.R
\name{boot}
\alias{boot}
\alias{boot.econet}
\title{boot: Bootstrap residuals with cross-sectional dependence}
\usage{
\method{boot}{econet}(
  object,
  hypothesis = c("lim", "het", "het_l", "het_r", "par", "par_split_with",
    "par_split_btw", "par_split_with_btw"),
  group = NULL,
  niter,
  weights = FALSE,
  delta = NULL,
  na.rm = FALSE,
  parallel = FALSE,
  cl,
  ...
)
}
\arguments{
\item{object}{an object of class \code{econet} obtained using 'NLLS' estimation.}

\item{hypothesis}{string. One of \code{c("lim","het", "het_l", "het_r", "par", "par_split_with", "par_split_btw", "par_split_with_btw")}.}

\item{group}{\code{NULL} or vector of positive integers specifying the indices for resampling within groups.}

\item{niter}{number of required iterations.}

\item{weights}{logical. It is \code{TRUE} if the object \code{object} was estimated using \code{weights}, and \code{FALSE} otherwise. Default \code{FALSE}.}

\item{delta}{Default is \code{NULL}.It has to be a number between zero (included) and one (excluded). When used, \code{econet} performs a constrained NLLS estimation. In this case, the estimated peer effect parameter, taken in absolute value, is forced to be higher than zero and lower than the spectral radius of \code{G}. Specifically, \code{delta} is a penalizing factor, decreasing the goodness of fit of the NLLS estimation, when the peer effect parameter approaches one of the two bounds. Observe that very high values of \code{delta} may cause NLLS estimation not to converge.}

\item{na.rm}{logical. Should missing values (including \code{NaN}) be removed?}

\item{parallel}{logical. It is \code{TRUE} if the user wants to make use of parallelization. Default is \code{FALSE}. Observe that this option is still in its beta version and it has been tested only for the Windows platform.}

\item{cl}{numeric. Number of cores to be used for parallelization.}

\item{...}{additional parameters}
}
\value{
a numeric vector containing bootstrapped standard errors  (see Anselin, 1990). If the procedure is not feasible, it returns a vector of NAs.
}
\description{
boot: Bootstrap residuals with cross-sectional dependence
}
\details{
For additional details, see the vignette (doi:10.18637/jss.v102.i08).
Warning: This function is available only when \code{net_dep} is run with \code{estimation == "NLLS"}
}
\examples{
\donttest{
# Load data
data("db_cosponsor")
data("G_alumni_111")
db_model_B <- db_cosponsor
G_model_B <- G_cosponsor_111
G_exclusion_restriction <- G_alumni_111
are_factors <- c("party", "gender", "nchair")
db_model_B[are_factors] <- lapply(db_model_B[are_factors], factor)

# Specify formula
f_model_B <- formula("les ~gender + party + nchair")

# Specify starting values
starting <- c(alpha = 0.23952,
              beta_gender1 = -0.22024,
              beta_party1 = 0.42947,
              beta_nchair1 = 3.09615,
              phi = 0.40038,
              unobservables = 0.07714)

# object Linear-in-means model
lim_model_B <- net_dep(formula = f_model_B, data = db_model_B,
                       G = G_model_B, model = "model_B", estimation = "NLLS",
                       hypothesis = "lim", endogeneity = TRUE, correction = "heckman",
                       first_step = "standard",
                       exclusion_restriction = G_exclusion_restriction,
                       start.val = starting)
# Bootstrap
# Warning: this may take a very long time to run.
# Decrease the number of iterations to reduce runtime.
# If you run econet on a Windows platform, you can try to set the
# argument parallel = TRUE. However note that this option is still
# in its beta version.
boot_lim_estimate <- boot(object = lim_model_B, hypothesis = "lim",
                          group = NULL, niter = 10, weights = FALSE)
boot_lim_estimate
}
}
\references{
Anselin, L., 1990, "Some robust approach to testing and estimation in spatial econometrics", Regional Science and Urban Economics, 20, 141-163.
}
\seealso{
\code{\link{net_dep}}
}
