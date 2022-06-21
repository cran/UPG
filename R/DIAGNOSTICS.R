#' @name UPG.Diag
#'
#' @title MCMC Diagnostics for \code{UPG.Probit}, \code{UPG.Logit}, \code{UPG.MNL} and \code{UPG.Binomial} objects using \code{coda}
#'
#' @description \code{UPG.Diag} computes a number of MCMC diagnostics based on the Markov chains that are contained in the model output returned by \code{UPG}.
#'
#' @param object an object of class \code{UPG.Probit}, \code{UPG.Logit}, \code{UPG.MNL} or \code{UPG.Binomial}.
#'
#' @return Returns a list containing effective sample size, effective sampling rate and inefficiency factors for each coefficient. In addition, maximum, minimum and median of these measures are returned.
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
#' # compute MCMC diagnostics
#' UPG.Diag(results.probit)
#'}
#'
#'@export
UPG.Diag           = function(object   = NULL    # estimated UPG object
                              ){


  if(!(class(object) %in% c("UPG.Probit","UPG.Logit","UPG.MNL","UPG.Binomial"))){

    stop("The input of UPG.Diag has to be an object generated using command UPG() in package UPG.")

  }


  if(object$inputs$model == "probit"){

    diagnostics = UPG.Diag.Probit(object)

  }

  if(object$inputs$model == "logit"){

    diagnostics = UPG.Diag.Logit(object)

  }

  if(object$inputs$model == "mnl"){

    diagnostics = UPG.Diag.MNL(object)

  }

  if(object$inputs$model == "binomial"){

    diagnostics = UPG.Diag.Binomial(object)

  }

  return(diagnostics)

}


#' @name UPG.Diag.Probit
#'
#' @title MCMC Diagnostics for UPG.Probit objects
#'
#' @description \code{UPG.Diag.Probit} computes inefficiency factors, effective sample size and effective sampling rate based on the posterior distributions in a \code{UPG.Probit} object.
#'
#' @param object an object of class \code{UPG.Probit}.
#'
#' @return Returns a list containing effective sample size, effective sampling rate and inefficiency factors for each coefficient.
#'
#' @author Gregor Zens
#'
UPG.Diag.Probit = function(object   = NULL    # estimated UPG object
                          ){

  # general info
  nsave = object$inputs$draws
  time  = object$runtime
  names = colnames(object$inputs$X)

  # ess, esr, ie
  ess   = coda::effectiveSize(coda::mcmc(object$posterior$beta))
  esr   = ess   / time
  ie    = nsave / ess

  names(ess) = names
  names(esr) = names
  names(ie)  = names

  # summaries ess
  ess.max = max(ess)
  ess.min = min(ess)
  ess.med = median(ess)

  # summaries esr
  esr.max = max(esr)
  esr.min = min(esr)
  esr.med = median(esr)

  # summaries ie
  ie.max  = max(ie)
  ie.min  = min(ie)
  ie.med  = median(ie)

  # create list
  diagnostics = list(

    inputs    = list(draws = nsave, time = time),

    summary   = list(

      ie      = list(max = ie.max,  min = ie.min,  median = ie.med),
      ess     = list(max = ess.max, min = ess.min, median = ess.med),
      esr     = list(max = esr.max, min = esr.min, median = esr.med)

    ),

    details   = list(ie  = ie, ess = ess, esr = esr)


  )


  return(diagnostics)

}

#' @name UPG.Diag.Logit
#'
#' @title MCMC Diagnostics for \code{UPG.Logit} objects
#'
#' @description \code{UPG.Diag.Logit} computes inefficiency factors, effective sample size and effective sampling rate based on the posterior distributions in a \code{UPG.Logit} object.
#'
#' @param object an object of class \code{UPG.Logit}.
#'
#' @return Returns a list containing effective sample size, effective sampling rate and inefficiency factors for each coefficient.
#'
#' @author Gregor Zens
#'
UPG.Diag.Logit = function(object   = NULL    # estimated UPG object
                          ){

  # general info
  nsave = object$inputs$draws
  time  = object$runtime
  names = colnames(object$inputs$X)

  # ess, esr, ie
  ess   = coda::effectiveSize(coda::mcmc(object$posterior$beta))
  esr   = ess   / time
  ie    = nsave / ess

  names(ess) = names
  names(esr) = names
  names(ie)  = names

  # summaries ess
  ess.max = max(ess)
  ess.min = min(ess)
  ess.med = median(ess)

  # summaries esr
  esr.max = max(esr)
  esr.min = min(esr)
  esr.med = median(esr)

  # summaries ie
  ie.max  = max(ie)
  ie.min  = min(ie)
  ie.med  = median(ie)

  # create list
  diagnostics = list(

    inputs    = list(draws = nsave, time = time),

    summary   = list(

      ie      = list(max = ie.max,  min = ie.min,  median = ie.med),
      ess     = list(max = ess.max, min = ess.min, median = ess.med),
      esr     = list(max = esr.max, min = esr.min, median = esr.med)

    ),

    details   = list(ie  = ie, ess = ess, esr = esr)


  )


  return(diagnostics)

}



#' @name UPG.Diag.Binomial
#'
#' @title MCMC Diagnostics for \code{UPG.Binomial} objects
#'
#' @description \code{UPG.Diag.Binomial} computes inefficiency factors, effective sample size and effective sampling rate based on the posterior distributions in a \code{UPG.Binomial} object.
#'
#' @param object an object of class \code{UPG.Binomial}.
#'
#' @return Returns a list containing effective sample size, effective sampling rate and inefficiency factors for each coefficient.
#'
#' @author Gregor Zens
#'
UPG.Diag.Binomial = function(object   = NULL    # estimated UPG object
                             ){

  # general info
  nsave = object$inputs$draws
  time  = object$runtime
  names = colnames(object$inputs$X)

  # ess, esr, ie
  ess   = coda::effectiveSize(coda::mcmc(object$posterior$beta))
  esr   = ess   / time
  ie    = nsave / ess

  names(ess) = names
  names(esr) = names
  names(ie)  = names

  # summaries ess
  ess.max = max(ess)
  ess.min = min(ess)
  ess.med = median(ess)

  # summaries esr
  esr.max = max(esr)
  esr.min = min(esr)
  esr.med = median(esr)

  # summaries ie
  ie.max  = max(ie)
  ie.min  = min(ie)
  ie.med  = median(ie)

  # create list
  diagnostics = list(

    inputs    = list(draws = nsave, time = time),

    summary   = list(

      ie      = list(max = ie.max,  min = ie.min,  median = ie.med),
      ess     = list(max = ess.max, min = ess.min, median = ess.med),
      esr     = list(max = esr.max, min = esr.min, median = esr.med)

    ),

    details   = list(ie  = ie, ess = ess, esr = esr)


  )


  return(diagnostics)

}


#' @name UPG.Diag.MNL
#'
#' @title MCMC Diagnostics for \code{UPG.MNL} objects
#'
#' @description \code{UPG.Diag.MNL} computes inefficiency factors, effective sample size and effective sampling rate based on the posterior distributions in a \code{UPG.MNL} object.
#'
#' @param object an object of class \code{UPG.MNL}.
#'
#' @return Returns a list containing effective sample size, effective sampling rate and inefficiency factors for each coefficient.
#'
#' @author Gregor Zens
#'
UPG.Diag.MNL = function(object   = NULL    # estimated UPG object
                        ){

  # general info
  nsave = object$inputs$draws
  time  = object$runtime
  names = colnames(object$inputs$X)
  K     = length(unique(object$inputs$y))

  # ess, esr, ie
  ess   = sapply(1:(K-1), function(kk)
                 coda::effectiveSize(coda::mcmc(object$posterior$beta[,,kk])))

  esr   = ess   / time
  ie    = nsave / ess

  rownames(ess) = names
  rownames(esr) = names
  rownames(ie)  = names

  # summaries ess
  ess.max = max(ess)
  ess.min = min(ess)
  ess.med = median(ess)

  # summaries esr
  esr.max = max(esr)
  esr.min = min(esr)
  esr.med = median(esr)

  # summaries ie
  ie.max  = max(ie)
  ie.min  = min(ie)
  ie.med  = median(ie)

  # create list
  diagnostics = list(

    inputs    = list(draws = nsave, time = time),

    summary   = list(

      ie      = list(max = ie.max,  min = ie.min,  median = ie.med),
      ess     = list(max = ess.max, min = ess.min, median = ess.med),
      esr     = list(max = esr.max, min = esr.min, median = esr.med)

    ),

    details   = list(ie  = ie, ess = ess, esr = esr),
    groups    = object$posterior$groups



  )


  return(diagnostics)

}
