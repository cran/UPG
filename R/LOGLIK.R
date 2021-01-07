#' @name logLik.UPG.Probit
#'
#' @title Compute log-likelihoods from UPG.Probit objects
#'
#' @description \code{logLik} can be used to compute log-likelihoods from \code{UPG.Probit} objects. The log-likelihood is based on the posterior mean of the coefficients and can be used for model selection when combined with,e.g., \code{BIC} or other model selection criteria.
#'
#' @param object an object of class \code{UPG.Probit}.
#' @param ... other logLik parameters.
#'
#' @return Returns a numeric of class \code{logLik} with attributes containing the number of estimated parameters and the number of observations.
#'
#' @seealso
#' \code{\link{summary.UPG.Probit}} to summarize the estimates of a discrete choice model from an \code{UPG.Probit} object and create tables.
#' \code{\link{plot.UPG.Probit}} to plot the results of a discrete choice model from an \code{UPG.Probit} object.
#' \code{\link{coef.UPG.Probit}} to extract coefficients from an \code{UPG.Probit} object.
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
#' # extract log-likelihood
#' ll.probit = logLik(results.probit)
#'
#' # compute BIC
#' BIC(ll.probit)
#'}
#' @method  logLik UPG.Probit
#'
#'@export
logLik.UPG.Probit =  function(object   = NULL,    # estimated UPG object
                              ...){

  means = apply(object$posterior$beta.post, 2, mean)
  fit   = object$inputs$X %*% means
  y     = object$inputs$y
  N     = nrow(object$inputs$X)
  K     = ncol(object$inputs$X)

  logl  = sum(y * log(pnorm(fit)) + (1-y) * log(1 - pnorm(fit)))

  attributes(logl) = list(nall = N, nobs = N, df = K)
  class(logl)      =  "logLik"

  return(logl)

}


#' @name logLik.UPG.Logit
#'
#' @title Compute log-likelihoods from UPG.Logit objects
#'
#' @description \code{logLik} can be used to compute log-likelihoods from \code{UPG.Logit} objects. The log-likelihood is based on the posterior mean of the coefficients and can be used for model selection when combined with,e.g., \code{BIC} or other model selection criteria.
#'
#' @param object an object of class \code{UPG.Logit}.
#' @param ... other logLik parameters.
#'
#' @return Returns a numeric of class \code{logLik} with attributes containing the number of estimated parameters and the number of observations.
#'
#' @seealso
#' \code{\link{summary.UPG.Logit}} to summarize the estimates of a discrete choice model from an \code{UPG.Logit} object and create tables.
#' \code{\link{plot.UPG.Logit}} to plot the results of a discrete choice model from an \code{UPG.Logit} object.
#' \code{\link{coef.UPG.Logit}} to extract coefficients from an \code{UPG.Logit} object.
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
#' # extract log-likelihood
#' ll.logit = logLik(results.logit)
#'
#' # compute BIC
#' BIC(ll.logit)
#'}
#' @method  logLik UPG.Logit
#'
#'@export
logLik.UPG.Logit  =  function(object   = NULL,    # estimated UPG object
                              ...){

  means = apply(object$posterior$beta.post, 2, mean)
  fit   = object$inputs$X %*% means
  y     = object$inputs$y
  N     = nrow(object$inputs$X)
  K     = ncol(object$inputs$X)

  logl  = sum(y * log(exp(fit) / (1 + exp(fit))) + (1-y) * log(1/(1+exp(fit))))

  attributes(logl) = list(nall = N, nobs = N, df = K)
  class(logl)      =  "logLik"

  return(logl)

}


#' @name logLik.UPG.MNL
#'
#' @title Compute log-likelihoods from UPG.MNL objects
#'
#' @description \code{logLik} can be used to compute log-likelihoods from \code{UPG.MNL} objects. The log-likelihood is based on the posterior mean of the coefficients and can be used for model selection when combined with,e.g., \code{BIC} or other model selection criteria.
#'
#' @param object an object of class \code{UPG.MNL}.
#' @param ... other logLik parameters.
#'
#' @return Returns a numeric of class \code{logLik} with attributes containing the number of estimated parameters and the number of observations.
#'
#' @seealso
#' \code{\link{summary.UPG.MNL}} to summarize the estimates of a discrete choice model from an \code{UPG.MNL} object and create tables.
#' \code{\link{plot.UPG.MNL}} to plot the results of a discrete choice model from an \code{UPG.MNL} object.
#' \code{\link{coef.UPG.MNL}} to extract coefficients from an \code{UPG.MNL} object.
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
#' # extract log-likelihood
#' ll.mnl = logLik(results.mnl)
#'
#' # compute BIC
#' BIC(ll.mnl)
#'}
#' @method  logLik UPG.MNL
#'
#'@export
logLik.UPG.MNL       =  function(object   = NULL,    # estimated UPG object
                                 ...){

  means = apply(object$posterior$beta.post, c(2,3), mean)
  fit   = object$inputs$X %*% means
  y     = object$inputs$y
  N     = nrow(object$inputs$X)
  noK   = length(unique(y)) # number of choices
  K     = ncol(object$inputs$X) * (noK-1) # number of parameters

  den   = log(rowSums(exp(fit)))

  logl  = sum(fit[y] - den)

  attributes(logl) = list(nall = N, nobs = N, df = K)
  class(logl)      =  "logLik"

  return(logl)

}


#' @name logLik.UPG.Binomial
#'
#' @title Compute log-likelihoods from UPG.Binomial objects
#'
#' @description \code{logLik} can be used to compute log-likelihoods from \code{UPG.Binomial} objects. The log-likelihood is based on the posterior mean of the coefficients and can be used for model selection when combined with,e.g., \code{BIC} or other model selection criteria.
#'
#' @param object an object of class \code{UPG.Binomial}.
#' @param ... other logLik parameters.
#'
#' @return Returns a numeric of class \code{logLik} with attributes containing the number of estimated parameters and the number of observations. Note that the number of observations in binomial models is equal to Ni and not equal to the sample size.
#'
#' @seealso
#' \code{\link{summary.UPG.Binomial}} to summarize the estimates of a discrete choice model from an \code{UPG.Binomial} object and create tables.
#' \code{\link{plot.UPG.Binomial}} to plot the results of a discrete choice model from an \code{UPG.Binomial} object.
#' \code{\link{coef.UPG.Binomial}} to extract coefficients from an \code{UPG.Binomial} object.
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
#' # extract log-likelihood
#' ll.binomial = logLik(results.binomial)
#'
#' # compute BIC
#' BIC(ll.binomial)
#'}
#' @method  logLik UPG.Binomial
#'
#'@export
logLik.UPG.Binomial  =  function(object   = NULL,    # estimated UPG object
                                 ...){

  means = apply(object$posterior$beta.post, 2, mean)
  fit   = object$inputs$X %*% means
  y     = object$inputs$y
  N     = sum(object$inputs$Ni) #note: don't use sample size following raftery 1995: BMS in Social Res.
  K     = ncol(object$inputs$X)
  Ni    = object$inputs$Ni

  logl  = sum(log(choose(Ni, y)) + y * log(exp(fit) / (1 + exp(fit))) + (Ni-y) * log(1/(1+exp(fit))))

  attributes(logl) = list(nall = N, nobs = N, df = K)
  class(logl)      =  "logLik"

  return(logl)

}
