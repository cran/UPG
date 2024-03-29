% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PLOT.R
\name{plot.UPG.Logit}
\alias{plot.UPG.Logit}
\title{Coefficient plots for UPG.Logit objects}
\usage{
\method{plot}{UPG.Logit}(
  x = NULL,
  ...,
  sort = FALSE,
  names = NULL,
  xlab = NULL,
  ylab = NULL,
  q = c(0.025, 0.975),
  include = NULL
)
}
\arguments{
\item{x}{an object of class \code{UPG.Logit}.}

\item{...}{other plot parameters.}

\item{sort}{a logical variable indicating whether the plotted coefficients should be sorted according to effect sizes. Default is FALSE.}

\item{names}{a character vector indicating names for the variables used in the plots.}

\item{xlab}{a character vector of length 1 indicating a title for the x-axis.}

\item{ylab}{a character vector of length 1 indicating a title for the y-axis.}

\item{q}{a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.}

\item{include}{can be used to plot only a subset of variables. Specify the columns of X that should be kept in the plot. See examples for further information.}
}
\value{
Returns a ggplot2 object.
}
\description{
\code{plot} generates plots from \code{UPG.Logit} objects using \code{ggplot2}. Coefficient plots show point estimates for all coefficients as well as their credible intervals.
}
\examples{
\donttest{
# estimate a logit model using example data
library(UPG)
data(lfp)
y = lfp[,1]
X = lfp[,-1]
results.logit = UPG(y = y, X = X, model = "logit")

# plot the results and sort coefficients by effect size
plot(results.logit, sort = TRUE)

# plot only variables 1 and 3 with custom names, credible intervals and axis labels
plot(results.logit,
     include  = c(1,3),
     names    = c("Custom 1", "Custom 2"),
     q        = c(0.1, 0.9),
     xlab     = c("Custom X"),
     ylab     = c("Custom Y"))
}
}
\seealso{
\code{\link{summary.UPG.Logit}} to summarize a \code{UPG.Logit} object and create tables.
\code{\link{predict.UPG.Logit}} to predict probabilities using a \code{UPG.Logit} object.
\code{\link{coef.UPG.Logit}} to extract coefficients from a \code{UPG.Logit} object.
}
\author{
Gregor Zens
}
