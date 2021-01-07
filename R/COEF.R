#' @name coef.UPG.Probit
#'
#' @title Extract coefficients from UPG.Probit objects
#'
#' @description \code{coef} can be used to extract posterior means and credible intervals based on posterior quantiles from \code{UPG.Probit} objects.
#'
#' @param object an object of class \code{UPG.Probit}.
#' @param ... other coef parameters.
#' @param q a numerical vector of length two holding the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#'
#' @return Returns a matrix containing posterior means and the desired credible interval.
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
#'
#' # extract posterior means and credible interval based on 0.025 and 0.975 quantiles
#' coef(results.probit, q = c(0.025, 0.975))
#'}
#' @method  coef UPG.Probit
#'
#'@export
coef.UPG.Probit = function(object,
                           ...,
                           q = c(0.025, 0.975)){


means = apply(object$posterior$beta.post, 2, mean)
lower = apply(object$posterior$beta.post, 2, quantile, q[1])
upper = apply(object$posterior$beta.post, 2, quantile, q[2])

coefs = cbind(lower, means, upper)

rownames(coefs) = colnames(object$inputs$X)
colnames(coefs) = c(paste0("Q",q[1] * 100), "Posterior Mean", paste0("Q",q[2] * 100))

return(coefs)


}


#' @name coef.UPG.Logit
#'
#' @title Extract coefficients from UPG.Logit objects
#'
#' @description \code{coef} can be used to extract posterior means and credible intervals based on posterior quantiles from \code{UPG.Logit} objects.
#'
#' @param object an object of class \code{UPG.Logit}.
#' @param ... other coef parameters.
#' @param q a numerical vector of length two holding the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#'
#' @return Returns a matrix containing posterior means and the desired credible interval.
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
#'
#' # extract posterior means and credible interval based on 0.025 and 0.975 quantiles
#' coef(results.logit, q = c(0.025, 0.975))
#'}
#' @method  coef UPG.Logit
#'
#'@export
coef.UPG.Logit  = function(object,
                           ...,
                           q = c(0.025, 0.975)
                           ){


  means = apply(object$posterior$beta.post, 2, mean)
  lower = apply(object$posterior$beta.post, 2, quantile, q[1])
  upper = apply(object$posterior$beta.post, 2, quantile, q[2])

  coefs = cbind(lower, means, upper)

  rownames(coefs) = colnames(object$inputs$X)
  colnames(coefs) = c(paste0("Q",q[1] * 100), "Posterior Mean", paste0("Q",q[2] * 100))

  return(coefs)


}



#' @name coef.UPG.MNL
#'
#' @title Extract coefficients from UPG.MNL objects
#'
#' @description \code{coef} can be used to extract posterior means and credible intervals based on posterior quantiles from \code{UPG.MNL} objects.
#'
#' @param object an object of class \code{UPG.MNL}.
#' @param ... other coef parameters.
#' @param q a numerical vector of length two holding the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#'
#' @return Returns a list containing posterior means and the desired credible interval.
#'
#' @seealso
#' \code{\link{summary.UPG.MNL}} to summarize the estimates of a discrete choice model from an \code{UPG.MNL} object and create tables.
#' \code{\link{predict.UPG.MNL}} to predict probabilities from a discrete choice model from an \code{UPG.MNL} object.
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
#'
#' # extract posterior means and credible interval based on 0.025 and 0.975 quantiles
#' coef(results.mnl, q = c(0.025, 0.975))
#'}
#' @method  coef UPG.MNL
#'
#'@export
coef.UPG.MNL    = function(object,
                           ...,
                           q = c(0.025, 0.975)
                           ){


  means = apply(object$posterior$beta.post, c(2,3), mean)
  lower = apply(object$posterior$beta.post, c(2,3), quantile, q[1])
  upper = apply(object$posterior$beta.post, c(2,3), quantile, q[2])

  rownames(means) = rownames(lower) = rownames(upper) = colnames(object$inputs$X)

  coefs = list(lower, means, upper, groups = object$posterior$groups)

  names(coefs)    = c(paste0("Q",q[1] * 100), "Posterior Mean", paste0("Q",q[2] * 100))

  return(coefs)


}



#' @name coef.UPG.Binomial
#'
#' @title Extract coefficients from UPG.Binomial objects
#'
#' @description \code{coef} can be used to extract posterior means and credible intervals based on posterior quantiles from \code{UPG.Binomial} objects.
#'
#' @param object an object of class \code{UPG.Binomial}.
#' @param ... other coef parameters.
#' @param q a numerical vector of length two holding the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#'
#' @return Returns a matrix containing posterior means and the desired credible interval.
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
#'
#' # extract posterior means and credible interval based on 0.025 and 0.975 quantiles
#' coef(results.binomial, q = c(0.025, 0.975))
#'}
#' @method  coef UPG.Binomial
#'
#'@export
coef.UPG.Binomial  = function(object,
                              ...,
                           q = c(0.025, 0.975)
                              ){


  means = apply(object$posterior$beta.post, 2, mean)
  lower = apply(object$posterior$beta.post, 2, quantile, q[1])
  upper = apply(object$posterior$beta.post, 2, quantile, q[2])

  coefs = cbind(lower, means, upper)

  rownames(coefs) = colnames(object$inputs$X)
  colnames(coefs) = c(paste0("Q",q[1] * 100), "Posterior Mean", paste0("Q",q[2] * 100))

  return(coefs)


}
