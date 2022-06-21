#' @name UPG
#'
#' @title Efficient MCMC Samplers for Bayesian probit regression and various logistic regression models
#'
#' @description \code{UPG} estimates Bayesian regression models for binary or categorical outcomes using samplers based on marginal data augmentation.
#'
#' @usage
#' UPG(y,
#'     X,
#'     model,
#'     Ni          = NULL,
#'     baseline    = NULL,
#'     draws       = 1000,
#'     burnin      = 1000,
#'     A0          = 4,
#'     B0          = 4,
#'     d0          = 2.5,
#'     D0          = 1.5,
#'     G0          = 100,
#'     verbose     = TRUE,
#'     gamma.boost = TRUE,
#'     delta.boost = TRUE,
#'     beta.start  = NULL)
#'
#' @param y a binary vector for probit and logit models. A character, factor or numeric vector for multinomial logit models. A numerical vector of the number of successes for the binomial model.
#' @param X a matrix of explanatory variables including an intercept in the first column. Rows are individuals, columns are variables.
#' @param model indicates the model to be estimated. \code{'probit'} for the probit model, \code{'logit'} for the logit model, \code{'mnl'} for the multinomial logit model or \code{'binomial'} for the binomial logit model.
#' @param Ni a vector containing the number of trials when estimating a binomial logit model.
#' @param baseline a string that can be used to change the baseline category in MNL models. Default baseline is the most commonly observed category.
#' @param draws number of saved Gibbs sampler iterations. Default is 1000 for illustration purposes, you should use more when estimating a model (e.g. 10,000).
#' @param burnin number of burned Gibbs sampler iterations. Default is 1000 for illustration purposes, you should use more when estimating a model (e.g. 2,000).
#' @param A0 prior variance for the intercept, 4 is the default.
#' @param B0 prior variance for the coefficients, 4 is the default.
#' @param d0 prior shape for working parameter delta, 2.5 is the default.
#' @param D0 prior rate for working parameter delta,  1.5 is the default.
#' @param G0 prior variance for working parameter gamma, 100 is the default.
#' @param verbose logical variable indicating whether progress should be printed during estimation.
#' @param gamma.boost logical variable indicating whether location-based parameter expansion boosting should be used.
#' @param delta.boost logical variable indicating whether scale-based parameter expansion boosting should be used.
#' @param beta.start provides starting values for beta (e.g. for use within Gibbs sampler)
#'
#' @return Depending on the estimated model, \code{UPG()} returns a \code{UPG.Probit}, \code{UPG.Logit}, \code{UPG.MNL} or \code{UPG.Binomial} object.
#'
#' @seealso
#' \code{\link{summary.UPG.Probit}} to summarize a \code{UPG.Probit} object and to create tables.
#' \code{\link{predict.UPG.Logit}} to predict probabilities using a \code{UPG.Logit} object.
#' \code{\link{plot.UPG.MNL}} to plot a \code{UPG.MNL} object.
#'
#' @author Gregor Zens
#'
#' @examples
#'
#' # load package
#' library(UPG)
#'
#' # estimate a probit model using example data
#' # warning: use more burn-ins, burnin = 100 is just used for demonstration purposes
#' data(lfp)
#' y = lfp[,1]
#' X = lfp[,-1]
#' results.probit = UPG(y = y, X = X, model = "probit", burnin = 100)
#'
#' # estimate a logit model using example data
#' # warning: use more burn-ins, burnin = 100 is just used for demonstration purposes
#' data(lfp)
#' y = lfp[,1]
#' X = lfp[,-1]
#' results.logit = UPG(y = y, X = X, model = "logit", burnin = 100)
#'
#' # estimate a MNL model using example data
#' # warning: use more burn-ins, burnin = 100 is just used for demonstration purposes
#' data(program)
#' y = program[,1]
#' X = program[,-1]
#' results.mnl = UPG(y = y, X = X, model = "mnl", burnin = 100)
#'
#' # estimate a binomial logit model using example data
#' # warning: use more burn-ins, burnin = 100 is just used for demonstration purposes
#' data(titanic)
#' y  = titanic[,1]
#' Ni = titanic[,2]
#' X  = titanic[,-c(1,2)]
#' results.binomial = UPG(y = y, X = X, Ni = Ni, model = "binomial", burnin = 100)
#'
#' @import ggplot2
#' @import knitr
#' @import stats
#' @import pgdraw
#' @import mnormt
#' @import matrixStats
#' @import coda
#' @importFrom utils flush.console setTxtProgressBar txtProgressBar
#'
#'@export
UPG = function(y,                     # dependent variable
               X,                     # design matrix
               model,                 # which model to estimate
               Ni          = NULL,    # number of trials (binomial only)
               baseline    = NULL,    # baseline category for MNL models
               draws       = 1000,    # number of saved draws
               burnin      = 1000,    # number of burn-in iterations
               A0          = 4,       # prior variance on intercepts
               B0          = 4,       # prior variance on remaining coefficients
               d0          = 2.5,     # prior scale of delta
               D0          = 1.5,     # prior rate  of delta
               G0          = 100,     # prior variance on gamma
               verbose     = TRUE,    # show progress bar
               gamma.boost = TRUE,    # use location-based expansion
               delta.boost = TRUE,    # use scale-based expansion
               beta.start  = NULL     # provide starting values for beta (for use in Gibbs sampler)
){

  if(verbose)          cat("Checking data & inputs ... \n")

  if(missing(y))        stop("Please provide a vector of dependent data y.")
  if(missing(X))        stop("Please provide a design matrix X.")
  if(missing(model))    stop("Please provide the model you want to estimate using the 'model' argument.")

  if(!(model %in% c("probit","logit","mnl","binomial"))) stop("'model' must be either 'probit', 'logit', 'mnl' or 'binomial'.")

  y = as.matrix(y, ncol=1)
  X = as.matrix(X, nrow=nrow(y))

  nsave = draws
  nburn = burnin

  # capture duplicated (and empty colnames)
  if(any(duplicated(colnames(X)))){stop("Please provide a design matrix X with unique column names.")}

  if(!is.null(beta.start)){beta.start = as.matrix(beta.start)}

  if(ncol(y)!=1)                       stop("y must be a vector.")
  if(!is.numeric(X))                   stop("Please provide X as numeric matrix.")
  if(is.null(X))                       stop("A design matrix X must be provided.")
  if(any(is.na(y), is.na(X)))          stop("Data contains NA values.")
  if(nrow(y) != nrow(X))               stop("y and X do not have the same number of rows.")
  if(sum(X[,1]) != nrow(X))            stop("The first column of X must be an intercept.")

  if(is.null(nburn))                   stop("Please provide the number of burn-ins using 'burnin'.")
  if(is.null(nsave))                   stop("Please provide the number of draws to be saved using 'draws'.")
  if(is.null(A0))                      stop("Please provide a scalar for 'A0'.")
  if(is.null(d0))                      stop("Please provide a scalar for 'd0'.")
  if(is.null(D0))                      stop("Please provide a scalar for 'D0'.")
  if(is.null(G0))                      stop("Please provide a scalar for 'G0'.")

  if(!is.null(beta.start)){

    if(any(is.na(beta.start))){stop("Starting value for beta contains NAs.")}
    if(!(class(beta.start)[1] %in% c("numeric","matrix"))){stop("Please provide a numeric vector or numeric matrix as starting value for beta.")}

  }

  if(verbose) cat("Initializing Gibbs Sampler ... \n")

  if(model == "probit"){

    if(!is.numeric(y))                   stop("Please provide y as numeric vector.")

    if(!(identical(as.integer(y), as.integer(as.logical(y))))) stop("The probit model requires a binary vector y.")

    if(!is.null(beta.start)){
      if(ncol(beta.start) != 1) stop(paste0("Starting value for beta needs to be a matrix with 1 column and ", ncol(X), " rows."))
      if(nrow(beta.start) != ncol(X)) stop(paste0("Starting value for beta needs to be a matrix with 1 column and ", ncol(X), " rows."))
    }

    if(verbose) cat("Simulating from posterior distribution ... \n")
    start.time = Sys.time()
    results    = upg.probit(y = y,
                            X = X,
                            nsave       = nsave,
                            nburn       = nburn,
                            d0          = d0,
                            D0          = D0,
                            G0          = G0,
                            B0          = B0,
                            A0          = A0,
                            gamma.boost = gamma.boost,
                            delta.boost = delta.boost,
                            beta.start  = beta.start,
                            verbose     = verbose)
    end.time   = Sys.time()
    output = list(beta   = results$beta)

    }

  if(model == "logit"){

    if(!is.numeric(y))                   stop("Please provide y as numeric vector.")

    if(!(identical(as.integer(y), as.integer(as.logical(y))))) stop("The logit model requires a binary vector y.")

    if(!is.null(beta.start)){
      if(ncol(beta.start) != 1) stop(paste0("Starting value for beta needs to be a matrix with 1 column and ", ncol(X), " rows."))
      if(nrow(beta.start) != ncol(X)) stop(paste0("Starting value for beta needs to be a matrix with 1 column and ", ncol(X), " rows."))

    }

    if(verbose) cat("Simulating from posterior distribution ... \n")

    start.time = Sys.time()
    results = upg.logit( y = y,
                         X = X,
                         nsave = nsave,
                         nburn = nburn,
                         d0 = d0,
                         D0 = D0,
                         G0 = G0,
                         B0 = B0,
                         A0 = A0,
                         gamma.boost = gamma.boost,
                         delta.boost = delta.boost,
                         beta.start  = beta.start,
                         verbose = verbose)
    end.time = Sys.time()


    output = list(beta   = results$beta)

  }

  if(model == "mnl"){

    if(length(unique(y)) == 2) warning("You have provided a binary outcome vector for a multinomial logit model. Are you sure that this is what you want?")

    if(!is.null(beta.start)){
      if(ncol(beta.start) != length(unique(y))) stop(paste0("Starting value for beta needs to be a matrix with ", length(unique(y)), " columns and ", ncol(X), " rows."))
      if(nrow(beta.start) != ncol(X)) stop(paste0("Starting value for beta needs to be a matrix with ", length(unique(y)), " columns and ", ncol(X), " rows."))
    }

    if(!is.null(baseline)){
    if(!(baseline %in% unique(y))){stop("Provided baseline is not part of the provided outcome vector y.")}
    }


    # MNL Baseline management

    # first, transform to factor
    y = factor(y)

    if(is.null(baseline)){

    # if no baseline is provided, use most common category
    tt       = table(y)
    baseline = names(tt[which.max(tt)])

    }

    pos.bl   = which(levels(y) == baseline)

    new.lvls = c(levels(y)[-pos.bl], baseline)

    # reorder factor with baseline last
    y.mnl    = factor(y, levels = new.lvls)

    # get groups for output, summary, plots, ...
    groups = levels(y.mnl)

    # transform to numeric categorical vector
    y.mnl = as.numeric(y.mnl)
    y.mnl = as.matrix(y.mnl)
    y.mnl = model.matrix(~as.factor(y.mnl) - 1)



    if(verbose) cat("Simulating from posterior distribution ... \n")

    start.time = Sys.time()
    results = upg.mnl(   y.matrix = y.mnl,
                         X        = X,
                         nsave    = nsave,
                         nburn    = nburn,
                         A0       = A0,
                         d0       = d0,
                         D0       = D0,
                         G0       = G0,
                         B0       = B0,
                         gamma.boost = gamma.boost,
                         delta.boost = delta.boost,
                         verbose  = verbose,
                         beta.start = beta.start)
    end.time  = Sys.time()


    output = list(beta   = results$beta,
                  groups      = groups)

  }

  if(model == "binomial"){

    if(!is.numeric(y))                   stop("Please provide y as numeric vector.")

    if(length(unique(y)) == 2) warning("You have provided a binary outcome vector for a binomial logit model. Are you sure that this is what you want?")
    if(is.null(Ni)) stop("The binomial logit model requires a separate vector containing the number of trials for each outcome.")
    if(is.null(Ni)) stop("The binomial logit model requires a separate vector containing the number of trials for each outcome.")
    if(any(Ni<y)) stop("Your data contains at least one case where the number of successes exceeds the number of trials.")
    Ni = as.matrix(Ni,ncol=1)
    if(nrow(y) != nrow(Ni))               stop("y and Ni do not have the same number of rows.")
    if(nrow(X) != nrow(Ni))               stop("X and Ni do not have the same number of rows.")
    if(!is.null(beta.start)){
      if(ncol(beta.start) != 1) stop(paste0("Starting value for beta needs to be a matrix with 1 column and ", ncol(X), " rows."))
      if(nrow(beta.start) != ncol(X)) stop(paste0("Starting value for beta needs to be a matrix with 1 column and ", ncol(X), " rows."))
    }

    if(verbose) cat("Simulating from posterior distribution ... \n")


    start.time = Sys.time()
    results = upg.binomial(y      = y,
                         X        = X,
                         nsave    = nsave,
                         nburn    = nburn,
                         A0       = A0,
                         d0       = d0,
                         D0       = D0,
                         G0       = G0,
                         B0       = B0,
                         gamma.boost = gamma.boost,
                         delta.boost = delta.boost,
                         verbose  = verbose,
                         beta.start = beta.start,
                         Ni       = Ni)
    end.time = Sys.time()


    output = list(beta = results$beta)

  }

  if(verbose) cat("\nSampling succesful!\n")


  #create output file

  input = list(y           = y,
               X           = X,
               model       = model,
               Ni          = Ni,
               baseline    = baseline,
               draws       = draws,
               burnin      = burnin,
               A0          = A0,
               B0          = B0,
               d0          = d0,
               D0          = D0,
               G0          = G0,
               verbose     = verbose,
               gamma.boost = gamma.boost,
               delta.boost = delta.boost,
               beta.start  = beta.start)




  # measure runtime of sampler
  runtime = as.numeric(difftime(end.time, start.time, units = "secs"))
  timecat = ifelse(runtime > 300,
                   paste0(round(runtime / 60, 2), " minutes."),
                   paste0(round(runtime, 2), " seconds."))



  # create output object
  if(verbose) cat("Saving output ...\n")
  UPG.output = list(posterior         = output,
                    inputs            = input,
                    runtime           = runtime)

  #s3 class attribute

  if(model == "probit")   class(UPG.output) = "UPG.Probit"
  if(model == "logit")    class(UPG.output) = "UPG.Logit"
  if(model == "mnl")      class(UPG.output) = "UPG.MNL"
  if(model == "binomial") class(UPG.output) = "UPG.Binomial"

  if(verbose) cat("Finished! Posterior simulation took",timecat,"\n")
  return(UPG.output)

}
