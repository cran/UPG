#' @name summary.UPG.Probit
#'
#' @title Estimation result summary and tables for UPG.Probit objects
#'
#' @description \code{summary} generates a summary of estimation results for \code{UPG.Probit} objects. Point estimates, estimated standard deviation as well as credible intervals for each variable are tabulated. In addition, an indicator quickly shows whether the credible interval includes zero or not. Moreover, LaTeX, HTML and pandoc tables can be quickly generated via \code{knitr}.
#'
#' @param object an object of class \code{UPG.Probit}.
#' @param q a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#' @param names a character vector indicating names for the variables used in the output.
#' @param digits number of digits to be included in output. Last digit will be rounded using \code{round}.
#' @param include can be used to summarize and tabulate only a subset of variables. Specify the columns of X that should be kept in the plot. See examples for further information.
#' @param table can be used to return a LaTeX table (\code{'latex'}), a Word table (\code{'pandoc'}) and HTML tables (\code{'html'}) via \code{knitr}. Include package "booktabs" in LaTeX preamble for LaTeX tables.
#' @param cap character vector that can be used to specify the table caption.
#' @param ... other summary parameters.
#'
#' @return Returns a \code{knitr_kable} object containing the summary table.
#'
#' @seealso
#' \code{\link{plot.UPG.Probit}} to plot a \code{UPG.Probit} object.
#' \code{\link{predict.UPG.Probit}} to predict probabilities using a \code{UPG.Probit} object.
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
#' # basic summary of regression results
#' summary(results.probit)
#'
#' # generate a LaTeX table with subset of variables and custom names
#' summary(results.probit,
#'         include=c(1,3),
#'         names=c("V. kept 1", "V. kept 3"),
#'         table="latex")
#'}
#' @method  summary UPG.Probit
#'
#'@export
summary.UPG.Probit = function(object = NULL,            #estimation output
                              ...,
                              q      = c(0.025, 0.975), #quantiles for credible interval
                              names  = NULL,            #variable names
                              digits = 2,               #round to digits
                              include= NULL,            #which variables to include? default:all (specify columns)
                              table  = NULL,            #create html or latex or pandoc table
                              cap    = NULL             #table caption
                       ){


  if(is.null(include))include = 1:ncol(object$inputs$X)
  if(is.null(names))  names   = colnames(object$inputs$X[,include,drop=F])
  if(is.null(names))  names   = paste0("Variable",1:ncol(object$inputs$X[,include,drop=F]))

  if(length(names) != length(include)) stop("Number of provided variable names does not match number of included variables.")

  #get summary statistics
  c.point  = apply(object$posterior$beta[,include,drop=F], 2, mean)
  c.upper  = apply(object$posterior$beta[,include,drop=F], 2, quantile, q[2])
  c.lower  = apply(object$posterior$beta[,include,drop=F], 2, quantile, q[1])
  c.sd     = apply(object$posterior$beta[,include,drop=F], 2, sd)

  #check which coefficients are "significant"
  c.sig = ifelse((c.upper > 0 & c.lower > 0) | (c.upper < 0 & c.lower < 0), "*", "")

  #make everything into table
  coefs = format(round(cbind(c.point, c.sd, c.lower, c.upper), digits), nsmall=digits)
  coefs = cbind(coefs, c.sig)

  hpd.upper = paste0("Q",100 * q[2])
  hpd.lower = paste0("Q",100 * q[1])

  header = c("Mean", "SD", hpd.lower, hpd.upper, paste0(100*(q[2]-q[1]), "% CI excl. 0"))

  df = data.frame(coefs)
  rownames(df) = names
  colnames(df) = header

  if(is.null(table)){

  result = kable(df, align = c("r","r","r","r","c"), row.names = T)


  } else {

    if(table=="latex"){
  result = kable(df[,-ncol(df)], align = c("r","r","r","r"),row.names = T,
                format = "latex",
                booktabs = T,
                linesep = "",
                caption = cap)
    }

    if(table=="pandoc"){
  result = kable(df[,-ncol(df)], align = c("r","r","r","r"),row.names = T,
                  format = "pandoc",
                  caption = cap)
    }

    if(table=="html"){
  result = kable(df[,-ncol(df)], align = c("r","r","r","r"),row.names = T,
                  format = "html",
                  caption = cap)
    }

  }

  # compute runtime, use in seconds or in minutes when > 5min
  runtime = object$runtime
  timecat = ifelse(runtime > 300,
                   paste0(round(runtime / 60, 2), " minutes."),
                   paste0(round(runtime, 2), " seconds."))

  cat("\n\n--- Bayesian Probit Results --- \n\n")
  cat("N =", length(object$inputs$y),"\n")
  cat("Analysis based on", object$inputs$draws, "posterior draws after\nan initial burn-in period of", object$inputs$burnin, "iterations.\n")
  cat("MCMC sampling took a total of",timecat,"\n")


  return(result)

}

#' @name summary.UPG.Logit
#'
#' @title Estimation results and tables for UPG.Logit objects
#'
#' @description \code{summary} generates a summary of estimation results for \code{UPG.Logit} objects. Point estimates, estimated standard deviation as well as credible intervals for each variable are tabulated. In addition, an indicator quickly shows whether the credible interval includes zero or not. Moreover, LaTeX, HTML and pandoc tables can be quickly generated via \code{knitr}.
#'
#' @param object an object of class \code{UPG.Logit}.
#' @param q a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#' @param names a character vector indicating names for the variables used in the output.
#' @param digits number of digits to be included in output. Last digit will be rounded using \code{round}.
#' @param include can be used to summarize and tabulate only a subset of variables. Specify the columns of X that should be kept in the plot. See examples for further information.
#' @param table can be used to return a LaTeX table (\code{'latex'}), a Word table (\code{'pandoc'}) and HTML tables (\code{'html'}) via \code{knitr}. Include package "booktabs" in LaTeX preamble for LaTeX tables.
#' @param cap character vector that can be used to specify the table caption.
#' @param ... other summary parameters.
#'
#' @return Returns a \code{knitr_kable} object containing the summary table.
#'
#' @seealso
#' \code{\link{plot.UPG.Logit}} to plot a \code{UPG.Logit} object.
#' \code{\link{predict.UPG.Logit}} to predict probabilities using a \code{UPG.Logit} object.
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
#' # basic summary of regression results
#' summary(results.logit)
#'
#' # generate a LaTeX table with subset of variables and custom names
#' summary(results.logit,
#'         include=c(1,3),
#'         names=c("V. kept 1", "V. kept 3"),
#'         table="latex")
#'}
#' @method  summary UPG.Logit
#'
#'@export
summary.UPG.Logit  = function(object = NULL,            #estimation output
                              ...,
                              q      = c(0.025, 0.975), #quantiles for credible interval
                              names  = NULL,            #variable names
                              digits = 2,               #round to digits
                              include= NULL,            #which variables to include? default:all (specify columns)
                              table  = NULL,            #create html or latex or pandoc table
                              cap    = NULL             #table caption
){


  if(is.null(include))include = 1:ncol(object$inputs$X)
  if(is.null(names))  names   = colnames(object$inputs$X[,include,drop=F])
  if(is.null(names))  names   = paste0("Variable",1:ncol(object$inputs$X[,include,drop=F]))

  if(length(names) != length(include)) stop("Number of provided variable names does not match number of included variables.")

  #get summary statistics
  c.point  = apply(object$posterior$beta[,include,drop=F], 2, mean)
  c.upper  = apply(object$posterior$beta[,include,drop=F], 2, quantile, q[2])
  c.lower  = apply(object$posterior$beta[,include,drop=F], 2, quantile, q[1])
  c.sd     = apply(object$posterior$beta[,include,drop=F], 2, sd)

  #check which coefficients are "significant"
  c.sig = ifelse((c.upper > 0 & c.lower > 0) | (c.upper < 0 & c.lower < 0), "*", "")

  #make everything into table
  coefs = format(round(cbind(c.point, c.sd, c.lower, c.upper), digits), nsmall=digits)
  coefs = cbind(coefs, c.sig)

  hpd.upper = paste0("Q",100 * q[2])
  hpd.lower = paste0("Q",100 * q[1])

  header = c("Mean", "SD", hpd.lower, hpd.upper, paste0(100*(q[2]-q[1]), "% CI excl. 0"))

  df = data.frame(coefs)
  rownames(df) = names
  colnames(df) = header

  if(is.null(table)){

    result = kable(df, align = c("r","r","r","r","c"), row.names = T)


  } else {

    if(table=="latex"){
      result = kable(df[,-ncol(df)], align = c("r","r","r","r"),row.names = T,
                     format = "latex",
                     booktabs = T,
                     linesep = "",
                     caption = cap)
    }

    if(table=="pandoc"){
      result = kable(df[,-ncol(df)], align = c("r","r","r","r"),row.names = T,
                     format = "pandoc",
                     caption = cap)
    }

    if(table=="html"){
      result = kable(df[,-ncol(df)], align = c("r","r","r","r"),row.names = T,
                     format = "html",
                     caption = cap)
    }

  }

  # compute runtime, use in seconds or in minutes when > 5min
  runtime = object$runtime
  timecat = ifelse(runtime > 300,
                   paste0(round(runtime / 60, 2), " minutes."),
                   paste0(round(runtime, 2), " seconds."))

  cat("\n\n--- Bayesian Logit Results --- \n\n")
  cat("N =", length(object$inputs$y),"\n")
  cat("Analysis based on", object$inputs$draws, "posterior draws after\nan initial burn-in period of", object$inputs$burnin, "iterations.\n")
  cat("MCMC sampling took a total of",timecat,"\n")


  return(result)

}




#' @name summary.UPG.MNL
#'
#' @title Estimation results and tables for UPG.MNL objects
#'
#' @description \code{summary} generates a summary of estimation results for \code{UPG.MNL} objects. Point estimates, estimated standard deviation as well as credible intervals for each variable are tabulated. In addition, an indicator quickly shows whether the credible interval includes zero or not. Moreover, LaTeX, HTML and pandoc tables can be quickly generated via \code{knitr}.
#'
#' @param object an object of class \code{UPG.MNL}.
#' @param q a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#' @param names a character vector indicating names for the variables used in the output.
#' @param groups  a character vector indicating names for the groups, excluding the baseline. The group names must correspond to the ordering in the dependent variable used for estimation.
#' @param digits number of digits to be included in output. Last digit will be rounded using \code{round}.
#' @param include can be used to summarize and tabulate only a subset of variables. Specify the columns of X that should be kept in the plot. See examples for further information.
#' @param table can be used to return a LaTeX table (\code{'latex'}), a Word table (\code{'pandoc'}) and HTML tables (\code{'html'}) via \code{knitr}. Include package "booktabs" in LaTeX preamble for LaTeX tables.
#' @param cap character vector that can be used to specify the table caption.
#' @param ... other summary parameters.
#'
#' @return Returns a \code{knitr_kable} object containing the summary table.
#'
#' @seealso
#' \code{\link{plot.UPG.MNL}} to plot a \code{UPG.MNL} object.
#' \code{\link{predict.UPG.MNL}} to predict probabilities using a \code{UPG.MNL} object.
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
#' # basic summary of regression results
#' summary(results.mnl)
#'
#' # generate a LaTeX table with subset of variables and custom names
#' summary(results.mnl,
#'         include=c(1,3),
#'         groups=c("Alpha","Beta"),
#'         names=c("V. kept 1", "V. kept 3"),
#'         table="latex")
#'}
#' @method  summary UPG.MNL
#'
#'@export
summary.UPG.MNL       = function(object = NULL,            #estimation output
                                 ...,
                                 q      = c(0.025, 0.975), #quantiles for credible interval
                                 groups = NULL,            #group names
                                 names  = NULL,            #variable names
                                 digits = 2,               #round to digits
                                 include= NULL,            #which variables to include? default:all (specify columns)
                                 table  = NULL,            #create html or latex or pandoc table
                                 cap    = NULL             #table caption
){


  if(is.null(include))include = 1:ncol(object$inputs$X)
  if(is.null(names))  names   = colnames(object$inputs$X[,include,drop=F])
  if(is.null(names))  names   = paste0("Variable",1:ncol(object$inputs$X[,include,drop=F]))

  if(length(names) != length(include)) stop("Number of provided variable names does not match number of included variables.")

  #get summary statistics
  c.point  = apply(object$posterior$beta[,include,,drop=F], c(2,3), mean)
  c.upper  = apply(object$posterior$beta[,include,,drop=F], c(2,3), quantile, q[2])
  c.lower  = apply(object$posterior$beta[,include,,drop=F], c(2,3), quantile, q[1])
  c.sd     = apply(object$posterior$beta[,include,,drop=F], c(2,3), sd)

  #check which coefficients are "significant"
  c.sig = ifelse((c.upper > 0 & c.lower > 0) | (c.upper < 0 & c.lower < 0), "*", "")

  #group names
  if(is.null(groups)) groups = paste("Category", paste0("'",object$posterior$groups[-ncol(c.point)],"'"))
  if(length(groups) != ncol(c.point)-1) stop("Wrong number of group names supplied. Need K-1 names where K is the number of choices.")

  #lose baseline
  c.sd    = c.sd[,-ncol(c.sd),drop=F]
  c.point = c.point[,-ncol(c.point),drop=F]
  c.sig   = c.sig[,-ncol(c.sig),drop=F]
  c.upper = c.upper[,-ncol(c.upper),drop=F]
  c.lower = c.lower[,-ncol(c.lower),drop=F]

  #create vectors for output
  c.point.mnl = c()
  c.sd.mnl    = c()
  c.sig.mnl   = c()
  c.upper.mnl = c()
  c.lower.mnl = c()

  for(ii in 1:ncol(c.point)){

    point.temp    = c("", format(round(c.point[,ii],digits), nsmall=digits),"")
    c.point.mnl   = c(c.point.mnl, point.temp)

    sd.temp       = c("", format(round(c.sd[,ii],digits), nsmall=digits),"")
    c.sd.mnl      = c(c.sd.mnl, sd.temp)

    sig.temp      = c("", c.sig[,ii],"")
    c.sig.mnl     = c(c.sig.mnl, sig.temp)

    upper.temp    = c("", format(round(c.upper[,ii],digits), nsmall=digits),"")
    c.upper.mnl   = c(c.upper.mnl, upper.temp)

    lower.temp    = c("", format(round(c.lower[,ii],digits), nsmall=digits),"")
    c.lower.mnl   = c(c.lower.mnl, lower.temp)


  }


  #make everything into table
  coefs        = cbind(c.point.mnl,c.sd.mnl,c.lower.mnl,c.upper.mnl, c.sig.mnl)
  hpd.upper    = paste0("Q",100*q[2])
  hpd.lower    = paste0("Q",100*q[1])
  header       = c("Mean", "SD", hpd.lower, hpd.upper, paste0(100*(q[2]-q[1]), "% CI excl. 0"))
  df           = data.frame(coefs)
  colnames(df) = header
  df           = df[-nrow(df),] # lose last line, unnecessary
  df           = as.matrix(df)

  #create row names
  row.temp     = c()

  for(ii in 1:(length(groups))){
    temp       = c(groups[ii],names,"")
    row.temp   = c(row.temp,temp)
  }

  row.temp     = row.temp[-length(row.temp)]
  rownames(df) = row.temp


  if(is.null(table)){

    result = kable(df, align = c("r","r","r","r","c"), row.names = T)


  } else {

    if(table=="latex"){
      result = kable(df[,-ncol(df)], align = c("r","r","r","r"),row.names = T,
                     format = "latex",
                     booktabs = T,
                     linesep = "",
                     caption = cap)
    }

    if(table=="pandoc"){
      result = kable(df[,-ncol(df)], align = c("r","r","r","r"),row.names = T,
                     format = "pandoc",
                     caption = cap)
    }

    if(table=="html"){
      result = kable(df[,-ncol(df)], align = c("r","r","r","r"),row.names = T,
                     format = "html",
                     caption = cap)
    }

  }

  # compute runtime, use in seconds or in minutes when > 5min
  runtime = object$runtime
  timecat = ifelse(runtime > 300,
                   paste0(round(runtime / 60, 2), " minutes."),
                   paste0(round(runtime, 2), " seconds."))

  cat("\n\n--- Bayesian Multinomial Logit Results --- \n\n")
  cat("N =", length(object$inputs$y),"\n")
  cat("Analysis based on", object$inputs$draws, "posterior draws after\nan initial burn-in period of", object$inputs$burnin, "iterations.\n")
  cat("MCMC sampling took a total of",timecat,"\n")
  cat("\n")
  cat("Category", paste0("'",object$posterior$groups[length(object$posterior$groups)],"'"), "is the baseline category.\n")

  return(result)

}






#' @name summary.UPG.Binomial
#'
#' @title Estimation results and tables for UPG.Binomial objects
#'
#' @description \code{summary} generates a summary of estimation results for \code{UPG.Binomial} objects. Point estimates, estimated standard deviation as well as credible intervals for each variable are tabulated. In addition, an indicator quickly shows whether the credible interval includes zero or not. Moreover, LaTeX, HTML and pandoc tables can be quickly generated via \code{knitr}.
#'
#' @param object an object of class \code{UPG.Binomial}.
#' @param q a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#' @param names a character vector indicating names for the variables used in the output.
#' @param digits number of digits to be included in output. Last digit will be rounded using \code{round}.
#' @param include can be used to summarize and tabulate only a subset of variables. Specify the columns of X that should be kept in the plot. See examples for further information.
#' @param table can be used to return a LaTeX table (\code{'latex'}), a Word table (\code{'pandoc'}) and HTML tables (\code{'html'}) via \code{knitr}. Include package "booktabs" in LaTeX preamble for LaTeX tables.
#' @param cap character vector that can be used to specify the table caption.
#' @param ... other summary parameters.
#'
#' @return Returns a \code{knitr_kable} object containing the summary table.
#'
#' @seealso
#' \code{\link{plot.UPG.Binomial}} to plot a \code{UPG.Binomial} object.
#' \code{\link{predict.UPG.Binomial}} to predict probabilities using a \code{UPG.Binomial} object.
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
#' # basic summary of regression results
#' summary(results.binomial)
#'
#' # generate a LaTeX table with subset of variables and custom names
#' summary(results.binomial,
#'         include=c(1,3),
#'         names=c("V. kept 1", "V. kept 3"),
#'         table="latex")
#'}
#' @method  summary UPG.Binomial
#'
#' @export
summary.UPG.Binomial  = function(object = NULL,            #estimation output
                                 ...,
                                 q      = c(0.025, 0.975), #quantiles for credible interval
                                 names  = NULL,            #variable names
                                 digits = 2,               #round to digits
                                 include= NULL,            #which variables to include? default:all (specify columns)
                                 table  = NULL,            #create html or latex or pandoc table
                                 cap    = NULL             #table caption
){


  if(is.null(include))include = 1:ncol(object$inputs$X)
  if(is.null(names))  names   = colnames(object$inputs$X[,include,drop=F])
  if(is.null(names))  names   = paste0("Variable",1:ncol(object$inputs$X[,include,drop=F]))

  if(length(names) != length(include)) stop("Number of provided variable names does not match number of included variables.")

  #get summary statistics
  c.point  = apply(object$posterior$beta[,include,drop=F], 2, mean)
  c.upper  = apply(object$posterior$beta[,include,drop=F], 2, quantile, q[2])
  c.lower  = apply(object$posterior$beta[,include,drop=F], 2, quantile, q[1])
  c.sd     = apply(object$posterior$beta[,include,drop=F], 2, sd)

  #check which coefficients are "significant"
  c.sig = ifelse((c.upper > 0 & c.lower > 0) | (c.upper < 0 & c.lower < 0), "*", "")

  #make everything into table
  coefs = format(round(cbind(c.point, c.sd, c.lower, c.upper), digits), nsmall=digits)
  coefs = cbind(coefs, c.sig)

  hpd.upper = paste0("Q",100 * q[2])
  hpd.lower = paste0("Q",100 * q[1])

  header = c("Mean", "SD", hpd.lower, hpd.upper, paste0(100*(q[2]-q[1]), "% CI excl. 0"))

  df = data.frame(coefs)
  rownames(df) = names
  colnames(df) = header

  if(is.null(table)){

    result = kable(df, align = c("r","r","r","r","c"), row.names = T)


  } else {

    if(table=="latex"){
      result = kable(df[,-ncol(df)], align = c("r","r","r","r"),row.names = T,
                     format = "latex",
                     booktabs = T,
                     linesep = "",
                     caption = cap)
    }

    if(table=="pandoc"){
      result = kable(df[,-ncol(df)], align = c("r","r","r","r"),row.names = T,
                     format = "pandoc",
                     caption = cap)
    }

    if(table=="html"){
      result = kable(df[,-ncol(df)], align = c("r","r","r","r"),row.names = T,
                     format = "html",
                     caption = cap)
    }

  }

  # compute runtime, use in seconds or in minutes when > 5min
  runtime = object$runtime
  timecat = ifelse(runtime > 300,
                   paste0(round(runtime / 60, 2), " minutes."),
                   paste0(round(runtime, 2), " seconds."))

  cat("\n\n--- Bayesian Binomial Logit Results --- \n\n")
  cat("N =", length(object$inputs$y),"\n")
  cat("Analysis based on", object$inputs$draws, "posterior draws after\nan initial burn-in period of", object$inputs$burnin, "iterations.\n")
  cat("MCMC sampling took a total of",timecat,"\n")

  return(result)

}
