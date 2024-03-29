% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PREDICT.R
\name{predict.UPG.MNL}
\alias{predict.UPG.MNL}
\title{Predicted probabilities from UPG.MNL objects}
\usage{
\method{predict}{UPG.MNL}(object = NULL, ..., newdata = NULL, q = c(0.025, 0.975))
}
\arguments{
\item{object}{an object of class \code{UPG.MNL}.}

\item{...}{other predict parameters.}

\item{newdata}{a matrix or a \code{data.frame} containing new explanatory data. The number of columns and the variable ordering must be the same as in the explanatory data used for estimation to generate valid predictions. If no new data is provided, \code{predict} will return predicted probabilities for the data used for estimating the model.}

\item{q}{a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.}
}
\value{
Returns a list containing posterior means of predicted probabilities as well as the desired credible interval.
}
\description{
\code{predict} generates predicted probabilities from a \code{UPG.MNL} object. In addition, credible intervals for these probabilities are computed. Probabilities can be predicted from the data used for estimating the model or for a new data set with the same structure.
}
\examples{
\donttest{
# estimate a multinomial logit model using example data
library(UPG)
data(program)
y = program[,1]
X = program[,-1]
results.mnl = UPG(y = y, X = X, model = "mnl")

# extract predicted probabilities
predict(results.mnl)
}
}
\seealso{
\code{\link{summary.UPG.MNL}} to summarize a \code{UPG.MNL} object and create tables.
\code{\link{plot.UPG.MNL}} to plot a \code{UPG.MNL} object.
\code{\link{coef.UPG.MNL}} to extract coefficients from a \code{UPG.MNL} object.
}
\author{
Gregor Zens
}
