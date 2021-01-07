#' @name print.UPG.Probit
#'
#' @title Print information for UPG.Probit objects
#'
#' @description \code{print} provides some basic information about an \code{UPG.Probit} object.
#'
#' @param x an object of class \code{UPG.Probit}.
#' @param ... other print parameters.
#'
#' @seealso
#' \code{\link{summary.UPG.Probit}} to summarize the estimates of a discrete choice model from an \code{UPG.Probit} object and create tables.
#' \code{\link{predict.UPG.Probit}} to predict probabilities from a discrete choice model from an \code{UPG.Probit} object.
#' \code{\link{plot.UPG.Probit}} to plot the results of a discrete choice model from an \code{UPG.Probit} object.
#'
#' @author Gregor Zens
#'
#' @examples
#' \donttest{
#' # estimate a probit model using example data
#' library(UPG)
#' data(lfp)
#' y = lfp[,1]
#' X = lfp[,-1]
#' results.probit = UPG(y = y, X = X, type = "probit", verbose=TRUE)
#' print(results.probit)
#'}
#' @method  print UPG.Probit
#'
#'@export
print.UPG.Probit = function(x,...){

cat("This is an UPG.Probit object. It contains the model inputs and the posterior distribution of a Bayesian probit model. The corresponding Gibbs sampler used information from", nrow(x$inputs$y),"data points and was iterated",x$inputs$nsave+x$inputs$nburn, "times, which took", x$runtime, "CPU seconds.")

  invisible(x)


}


#' @name print.UPG.Logit
#'
#' @title Print information for UPG.Logit objects
#'
#' @description \code{print} provides some basic information about an \code{UPG.Logit} object.
#'
#' @param x an object of class \code{UPG.Logit}.
#' @param ... other print parameters.
#'
#' @seealso
#' \code{\link{summary.UPG.Logit}} to summarize the estimates of a discrete choice model from an \code{UPG.Logit} object and create tables.
#' \code{\link{predict.UPG.Logit}} to predict probabilities from a discrete choice model from an \code{UPG.Logit} object.
#' \code{\link{plot.UPG.Logit}} to plot the results of a discrete choice model from an \code{UPG.Logit} object.
#'
#' @author Gregor Zens
#'
#' @examples
#' \donttest{
#' # estimate a logit model using example data
#' library(UPG)
#' data(lfp)
#' y = lfp[,1]
#' X = lfp[,-1]
#' results.logit = UPG(y = y, X = X, type = "logit", verbose=TRUE)
#' print(results.logit)
#'}
#' @method  print UPG.Logit
#'
#'@export
print.UPG.Logit = function(x,...){

  cat("This is an UPG.Logit object. It contains the model inputs and the posterior distribution of a Bayesian logit model. The corresponding Gibbs sampler used information from", nrow(x$inputs$y),"data points and was iterated",x$inputs$nsave+x$inputs$nburn, "times, which took", x$runtime, "CPU seconds.")

  invisible(x)


}


#' @name print.UPG.MNL
#'
#' @title Print information for UPG.MNL objects
#'
#' @description \code{print} provides some basic information about an \code{UPG.MNL} object.
#'
#' @param x an object of class \code{UPG.MNL}.
#' @param ... other print parameters.
#'
#' @seealso
#' \code{\link{summary.UPG.MNL}} to summarize the estimates of a discrete choice model from an \code{UPG.MNL} object and create tables.
#' \code{\link{predict.UPG.MNL}} to predict probabilities from a discrete choice model from an \code{UPG.Logit} object.
#' \code{\link{plot.UPG.MNL}} to plot the results of a discrete choice model from an \code{UPG.MNL} object.
#'
#' @author Gregor Zens
#'
#' @examples
#' \donttest{
#' # estimate a multinomial logit model using example data
#' library(UPG)
#' data(program)
#' y = program[,1]
#' X = program[,-1]
#' results.mnl = UPG(y = y, X = X, type = "mnl")
#' print(results.mnl)
#'}
#' @method  print UPG.MNL
#'
#'@export
print.UPG.MNL = function(x,...){

  cat("This is an UPG.MNL object. It contains the model inputs and the posterior distribution of a Bayesian multinomial logit model. The corresponding Gibbs sampler used information from", nrow(x$inputs$y),"data points and was iterated",x$inputs$nsave+x$inputs$nburn, "times, which took", x$runtime, "CPU seconds.")

  invisible(x)


}


#' @name print.UPG.Binomial
#'
#' @title Print information for UPG.Binomial objects
#'
#' @description \code{print} provides some basic information about an \code{UPG.Binomial} object.
#'
#' @param x an object of class \code{UPG.Binomial}.
#' @param ... other print parameters.
#'
#' @seealso
#' \code{\link{summary.UPG.Binomial}} to summarize the estimates of a discrete choice model from an \code{UPG.Binomial} object and create tables.
#' \code{\link{predict.UPG.Binomial}} to predict probabilities from a discrete choice model from an \code{UPG.Binomial} object.
#' \code{\link{plot.UPG.Binomial}} to plot the results of a discrete choice model from an \code{UPG.Binomial} object.
#'
#' @author Gregor Zens
#'
#' @examples
#' \donttest{
#' # estimate a binomial logit model using example data
#' library(UPG)
#' data(titanic)
#' y  = titanic[,1]
#' Ni = titanic[,2]
#' X  = titanic[,-c(1,2)]
#' results.binomial = UPG(y = y, X = X, Ni = Ni, type = "binomial")
#' print(results.binomial)
#'}
#' @method  print UPG.Binomial
#'
#'@export
print.UPG.Binomial = function(x,...){

  cat("This is an UPG.Binomial object. It contains the model inputs and the posterior distribution of a Bayesian binomial logit model. The corresponding Gibbs sampler used information from", nrow(x$inputs$y),"data points and was iterated",x$inputs$nsave+x$inputs$nburn, "times, which took", x$runtime, "CPU seconds.")


  invisible(x)

}
