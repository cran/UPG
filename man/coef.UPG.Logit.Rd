% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/COEF.R
\name{coef.UPG.Logit}
\alias{coef.UPG.Logit}
\title{Extract coefficients from UPG.Logit objects}
\usage{
\method{coef}{UPG.Logit}(object, ..., q = c(0.025, 0.975))
}
\arguments{
\item{object}{an object of class \code{UPG.Logit}.}

\item{...}{other coef parameters.}

\item{q}{a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.}
}
\value{
Returns a matrix containing posterior means and the desired credible interval.
}
\description{
\code{coef} can be used to extract posterior means and credible intervals based on posterior quantiles from \code{UPG.Logit} objects.
}
\examples{
\donttest{
# estimate a logit model using example data
library(UPG)
data(lfp)
y = lfp[,1]
X = lfp[,-1]
results.logit = UPG(y = y, X = X, model = "logit")

# extract posterior means and credible interval based on 0.025 and 0.975 quantiles
coef(results.logit, q = c(0.025, 0.975))
}
}
\seealso{
\code{\link{summary.UPG.Logit}} to summarize a \code{UPG.Logit} object and create tables.
\code{\link{predict.UPG.Logit}} to predict probabilities using a \code{UPG.Logit} object.
\code{\link{plot.UPG.Logit}} to plot a \code{UPG.Logit} object.
}
\author{
Gregor Zens
}
