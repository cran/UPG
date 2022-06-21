#' @name predict.UPG.Probit
#'
#' @title Predicted probabilities from UPG.Probit objects
#'
#' @description \code{predict} generates predicted probabilities from a \code{UPG.Probit} object. In addition, credible intervals for these probabilities are computed. Probabilities can be predicted from the data used for estimating the model or for a new data set with the same structure.
#'
#' @param object an object of class \code{UPG.Probit}.
#' @param newdata a matrix or a \code{data.frame} containing new explanatory data. The number of columns and the variable ordering must be the same as in the explanatory data used for estimation to generate valid predictions. If no new data is provided, \code{predict} will return predicted probabilities for the data used for estimating the model.
#' @param q a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#' @param ... other predict parameters.
#'
#' @return Returns a list containing posterior means of predicted probabilities as well as the desired credible interval.
#'
#' @seealso
#' \code{\link{summary.UPG.Probit}} to summarize a \code{UPG.Probit} object and create tables.
#' \code{\link{plot.UPG.Probit}} to plot a \code{UPG.Probit} object.
#' \code{\link{coef.UPG.Probit}} to extract coefficients from a \code{UPG.Probit} object.
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
#' results.probit = UPG(y = y, X = X, model = "probit")
#'
#' # extract predicted probabilities
#' predict(results.probit)
#'}
#' @method  predict UPG.Probit
#'
#'@export
predict.UPG.Probit = function(object   = NULL,    # estimated UPG object
                              ...,
                              newdata = NULL,     # new data to predict from (if NULL, estimation data is used)
                              q = c(0.025, 0.975) # quantiles used for credible intervals
                              ){

  if(is.null(newdata)){
    Z = as.matrix(object$inputs$X)
  } else {
    Z = as.matrix(newdata)
  }

  if(!is.null(newdata)){
  if(ncol(newdata) != ncol(object$inputs$X)) stop("Number of variables of prediction data does not match number of variables used in estimation.")
  }

  # shell for predicted probabilities
  ppred = array(NA, c(object$inputs$draws, nrow(Z)))

  # compute predicted probabilities
  for(ii in 1:object$inputs$draws){

    ppred[ii,] = pnorm(Z %*% object$posterior$beta[ii,])

  }

  means = apply(ppred, 2, mean)

  lower = apply(ppred, 2, quantile, q[1])

  upper = apply(ppred, 2, quantile, q[2])


  prediction = list(lower, means, upper)


  names(prediction) = c(paste0("Q", 100 * q[1]), "Posterior mean", paste0("Q", 100 * q[2]))

  return(prediction)

}




#' @name predict.UPG.Logit
#'
#' @title Predicted probabilities from UPG.Logit objects
#'
#' @description \code{predict} generates predicted probabilities from a \code{UPG.Logit} object. In addition, credible intervals for these probabilities are computed. Probabilities can be predicted from the data used for estimating the model or for a new data set with the same structure.
#'
#' @param object an object of class \code{UPG.Logit}.
#' @param newdata a matrix or a \code{data.frame} containing new explanatory data. The number of columns and the variable ordering must be the same as in the explanatory data used for estimation to generate valid predictions. If no new data is provided, \code{predict} will return predicted probabilities for the data used for estimating the model.
#' @param q a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#' @param ... other predict parameters.
#'
#' @return Returns a list containing posterior means of predicted probabilities as well as the desired credible interval.
#'
#' @seealso
#' \code{\link{summary.UPG.Logit}} to summarize a \code{UPG.Logit} object and create tables.
#' \code{\link{plot.UPG.Logit}} to plot a \code{UPG.Logit} object.
#' \code{\link{coef.UPG.Logit}} to extract coefficients from a \code{UPG.Logit} object.
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
#' results.logit = UPG(y = y, X = X, model = "logit")
#'
#' # extract predicted probabilities
#' predict(results.logit)
#'}
#' @method  predict UPG.Logit
#'
#'@export
predict.UPG.Logit = function(object   = NULL,    # estimated UPG object
                             ...,
                             newdata = NULL,     # new data to predict from (if NULL, estimation data is used)
                             q = c(0.025, 0.975) # quantiles used for credible intervals
){

  if(is.null(newdata)){
    Z = as.matrix(object$inputs$X)
  } else {
    Z = as.matrix(newdata)
  }

  if(!is.null(newdata)){
    if(ncol(newdata) != ncol(object$inputs$X)) stop("Number of variables of prediction data does not match number of variables used in estimation.")
  }


  # shell for predicted probabilities
  ppred = array(NA, c(object$inputs$draws, nrow(Z)))

  # compute predicted probabilities
  for(ii in 1:object$inputs$draws){

    ppred[ii,] = 1 / (1 + exp(-Z %*% object$posterior$beta[ii,]))

  }

  means = apply(ppred, 2, mean)

  lower = apply(ppred, 2, quantile, q[1])

  upper = apply(ppred, 2, quantile, q[2])


  prediction = list(lower, means, upper)


  names(prediction) = c(paste0("Q", 100 * q[1]), "Posterior mean", paste0("Q", 100 * q[2]))

  return(prediction)

}



#' @name predict.UPG.MNL
#'
#' @title Predicted probabilities from UPG.MNL objects
#'
#' @description \code{predict} generates predicted probabilities from a \code{UPG.MNL} object. In addition, credible intervals for these probabilities are computed. Probabilities can be predicted from the data used for estimating the model or for a new data set with the same structure.
#'
#' @param object an object of class \code{UPG.MNL}.
#' @param newdata a matrix or a \code{data.frame} containing new explanatory data. The number of columns and the variable ordering must be the same as in the explanatory data used for estimation to generate valid predictions. If no new data is provided, \code{predict} will return predicted probabilities for the data used for estimating the model.
#' @param q a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#' @param ... other predict parameters.
#'
#' @return Returns a list containing posterior means of predicted probabilities as well as the desired credible interval.
#'
#' @seealso
#' \code{\link{summary.UPG.MNL}} to summarize a \code{UPG.MNL} object and create tables.
#' \code{\link{plot.UPG.MNL}} to plot a \code{UPG.MNL} object.
#' \code{\link{coef.UPG.MNL}} to extract coefficients from a \code{UPG.MNL} object.
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
#' results.mnl = UPG(y = y, X = X, model = "mnl")
#'
#' # extract predicted probabilities
#' predict(results.mnl)
#'}
#' @method  predict UPG.MNL
#'
#'@export
predict.UPG.MNL   = function(object   = NULL,    # estimated UPG object
                             ...,
                             newdata = NULL,     # new data to predict from (if NULL, estimation data is used)
                             q = c(0.025, 0.975) # quantiles used for credible intervals
){

  if(is.null(newdata)){
    Z = as.matrix(object$inputs$X)
  } else {
    Z = as.matrix(newdata)
  }

  if(!is.null(newdata)){
    if(ncol(newdata) != ncol(object$inputs$X)) stop("Number of variables of prediction data does not match number of variables used in estimation.")
  }


  # prediction shell
  ppred = array(NA, c(object$inputs$draws, nrow(Z), length(unique(object$inputs$y))))

  # compute predictions
  for(ii in 1:object$inputs$draws){

    ppred[ii,,] = exp(Z %*% object$posterior$beta[ii,,])
    ppred[ii,,] = ppred[ii,,] / rowSums(ppred[ii,,])

  }

  means = apply(ppred, c(2,3), mean)

  lower = apply(ppred, c(2,3), quantile, q[1])

  upper = apply(ppred, c(2,3), quantile, q[2])

  prediction = list(lower, means, upper, groups = object$posterior$groups)


  names(prediction) = c(paste0("Q", 100 * q[1]), "Posterior mean", paste0("Q", 100 * q[2]), "Groups")

  return(prediction)

}













#' @name predict.UPG.Binomial
#'
#' @title Predicted probabilities from UPG.Binomial objects
#'
#' @description \code{predict} generates predicted probabilities from a \code{UPG.Binomial} object. In addition, credible intervals for these probabilities are computed. Probabilities can be predicted from the data used for estimating the model or for a new data set with the same structure.
#'
#' @param object an object of class \code{UPG.Binomial}.
#' @param newdata a matrix or a \code{data.frame} containing new explanatory data. The number of columns and the variable ordering must be the same as in the explanatory data used for estimation to generate valid predictions. If no new data is provided, \code{predict} will return predicted probabilities for the data used for estimating the model.
#' @param q a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#' @param ... other predict parameters.
#'
#' @return Returns a list containing posterior means of predicted probabilities as well as the desired credible interval.
#'
#' @seealso
#' \code{\link{summary.UPG.Binomial}} to summarize a \code{UPG.Binomial} object and create tables.
#' \code{\link{plot.UPG.Binomial}} to plot a \code{UPG.Binomial} object.
#' \code{\link{coef.UPG.Binomial}} to extract coefficients from a \code{UPG.Binomial} object.
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
#' results.binomial = UPG(y = y, X = X, Ni = Ni, model = "binomial")
#'
#' # extract predicted probabilities
#' predict(results.binomial)
#'}
#' @method  predict UPG.Binomial
#'
#'@export
predict.UPG.Binomial = function(object   = NULL,    # estimated UPG object
                                ...,
                                newdata = NULL,     # new data to predict from (if NULL, estimation data is used)
                                q = c(0.025, 0.975) # quantiles used for credible intervals
){

  if(is.null(newdata)){
    Z = as.matrix(object$inputs$X)
  } else {
    Z = as.matrix(newdata)
  }

  if(!is.null(newdata)){
    if(ncol(newdata) != ncol(object$inputs$X)) stop("Number of variables of prediction data does not match number of variables used in estimation.")
  }


  # shell for predicted probabilities
  ppred = array(NA, c(object$inputs$draws, nrow(Z)))

  # compute predicted probabilities
  for(ii in 1:object$inputs$draws){

    ppred[ii,] = 1 / (1 + exp(-Z %*% object$posterior$beta[ii,]))

  }

  means = apply(ppred, 2, mean)

  lower = apply(ppred, 2, quantile, q[1])

  upper = apply(ppred, 2, quantile, q[2])


  prediction = list(lower, means, upper)



  names(prediction) = c(paste0("Q", 100 * q[1]), "Posterior mean", paste0("Q", 100 * q[2]))

  return(prediction)

}
