% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PRINT.R
\name{print.UPG.MNL}
\alias{print.UPG.MNL}
\title{Print information for UPG.MNL objects}
\usage{
\method{print}{UPG.MNL}(x, ...)
}
\arguments{
\item{x}{an object of class \code{UPG.MNL}.}

\item{...}{other print parameters.}
}
\description{
\code{print} provides some basic information about a \code{UPG.MNL} object.
}
\examples{
\donttest{
# estimate a multinomial logit model using example data
library(UPG)
data(program)
y = program[,1]
X = program[,-1]
results.mnl = UPG(y = y, X = X, model = "mnl")
print(results.mnl)
}
}
\seealso{
\code{\link{summary.UPG.MNL}} to summarize a \code{UPG.MNL} object and create tables.
\code{\link{predict.UPG.MNL}} to predict probabilities using a \code{UPG.MNL} object.
\code{\link{plot.UPG.MNL}} to plot a \code{UPG.MNL} object.
}
\author{
Gregor Zens
}
