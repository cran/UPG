#' @name plot.UPG.Probit
#' @title Coefficient plots for UPG.Probit objects
#'
#' @description \code{plot} generates plots from \code{UPG.Probit} objects using \code{ggplot2}. Coefficient plots show point estimates for all coefficients as well as their credible intervals.
#'
#' @param x an object of class \code{UPG.Probit}.
#' @param sort a logical variable indicating whether the plotted coefficients should be sorted according to effect sizes. Default is FALSE.
#' @param names a character vector indicating names for the variables used in the plots.
#' @param xlab  a character vector of length 1 indicating a title for the x-axis.
#' @param ylab  a character vector of length 1 indicating a title for the y-axis.
#' @param q a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#' @param include can be used to plot only a subset of variables. Specify the columns of X that should be kept in the plot. See examples for further information.
#' @param ... other plot parameters.
#'
#' @return Returns a ggplot2 object.
#'
#' @seealso
#' \code{\link{summary.UPG.Probit}} to summarize a \code{UPG.Probit} object and create tables.
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
#' # plot the results and sort coefficients by effect size
#' plot(results.probit, sort = TRUE)
#'
#' # plot only variables 1 and 3 with custom names, credible intervals and axis labels
#' plot(results.probit,
#'      include  = c(1, 3),
#'      names    = c("Custom 1", "Custom 2"),
#'      q        = c(0.1, 0.9),
#'      xlab     = c("Custom X"),
#'      ylab     = c("Custom Y"))
#'}
#' @method  plot UPG.Probit
#'
#'@export
plot.UPG.Probit = function(x         = NULL,
                           ...,
                           sort      = FALSE,           # sort coefficients by average effect size
                           names     = NULL,            # provide names for variables, alternatively
                           xlab      = NULL,            # provide x axis label
                           ylab      = NULL,            # provide y axis label
                           q         = c(0.025, 0.975), # credible intervals
                           include   = NULL             # which variables to include? default:all (numeric vector)
                           ){



  if(is.null(include)) include = 1:ncol(x$inputs$X)
  if(is.null(names))   names = colnames(x$inputs$X[,include,drop=F])
  if(is.null(names))   names = paste0("Variable", 1:ncol(x$inputs$X[,include,drop=F]))
  if(length(names) != length(include)) stop("Number of provided variable names does not match number of included variables.")

  #create some global variables such that r cmd check is happy
  variable = iter = value = NULL

  c.point  = apply(x$posterior$beta[,include,drop=F], 2, mean)
  c.upper  = apply(x$posterior$beta[,include,drop=F], 2, quantile, q[2])
  c.lower  = apply(x$posterior$beta[,include,drop=F], 2, quantile, q[1])

  plot.df  = data.frame(c.point,c.upper,c.lower,names)

  if(nrow(plot.df) > 1 & sort){

    plot.df$names = factor(plot.df$names, levels = plot.df$names[order(plot.df$c.point)])

  }

  if(!sort){

    # plot in order of appearance in X
    plot.df$names = factor(plot.df$names, levels = rev(plot.df$names))

  }

  #axis labeling, be careful b/c of coord_flip
  if(is.null(ylab)) ylab = ""
  if(is.null(xlab)) xlab = "Posterior Estimate"

  final = ggplot(plot.df, aes(x=names, y=c.point)) +
           geom_errorbar(ymin=c.lower, ymax=c.upper, width=0, col="grey60") +
           geom_point() +
           geom_point(data = plot.df, aes(x=names, y=c.upper), size=-1) + #for correct scaling
           geom_point(data = plot.df, aes(x=names, y=c.lower), size=-1) + #for correct scaling
           theme_bw()   +
           geom_hline(yintercept = 0,lty=2) +
           coord_flip() +
           xlab(ylab)   +
           ylab(xlab)

  return(final)

}




#' @name plot.UPG.Logit
#' @title Coefficient plots for UPG.Logit objects
#'
#' @description \code{plot} generates plots from \code{UPG.Logit} objects using \code{ggplot2}. Coefficient plots show point estimates for all coefficients as well as their credible intervals.
#'
#' @param x an object of class \code{UPG.Logit}.
#' @param sort a logical variable indicating whether the plotted coefficients should be sorted according to effect sizes. Default is FALSE.
#' @param names a character vector indicating names for the variables used in the plots.
#' @param xlab  a character vector of length 1 indicating a title for the x-axis.
#' @param ylab  a character vector of length 1 indicating a title for the y-axis.
#' @param q a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#' @param include can be used to plot only a subset of variables. Specify the columns of X that should be kept in the plot. See examples for further information.
#' @param ... other plot parameters.
#'
#' @return Returns a ggplot2 object.
#'
#' @seealso
#' \code{\link{summary.UPG.Logit}} to summarize a \code{UPG.Logit} object and create tables.
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
#' # plot the results and sort coefficients by effect size
#' plot(results.logit, sort = TRUE)
#'
#' # plot only variables 1 and 3 with custom names, credible intervals and axis labels
#' plot(results.logit,
#'      include  = c(1,3),
#'      names    = c("Custom 1", "Custom 2"),
#'      q        = c(0.1, 0.9),
#'      xlab     = c("Custom X"),
#'      ylab     = c("Custom Y"))
#'}
#' @method  plot UPG.Logit
#'
#'@export
plot.UPG.Logit  = function(x         = NULL,
                           ...,
                           sort      = FALSE,           # sort coefficients by average effect size
                           names     = NULL,            # provide names for variables, alternatively
                           xlab      = NULL,            # provide x axis label
                           ylab      = NULL,            # provide y axis label
                           q         = c(0.025, 0.975), # credible intervals
                           include   = NULL             # which variables to include? default:all (numeric vector)
){



  if(is.null(include)) include = 1:ncol(x$inputs$X)
  if(is.null(names))   names = colnames(x$inputs$X[,include,drop=F])
  if(is.null(names))   names = paste0("Variable", 1:ncol(x$inputs$X[,include,drop=F]))
  if(length(names) != length(include)) stop("Number of provided variable names does not match number of included variables.")

  #create some global variables such that r cmd check is happy
  variable = iter = value = NULL

  c.point  = apply(x$posterior$beta[,include,drop=F], 2, mean)
  c.upper  = apply(x$posterior$beta[,include,drop=F], 2, quantile, q[2])
  c.lower  = apply(x$posterior$beta[,include,drop=F], 2, quantile, q[1])

  plot.df  = data.frame(c.point,c.upper,c.lower,names)

  if(nrow(plot.df) > 1 & sort){

    plot.df$names = factor(plot.df$names, levels = plot.df$names[order(plot.df$c.point)])

  }

  if(!sort){

    # plot in order of appearance in X
    plot.df$names = factor(plot.df$names, levels = rev(plot.df$names))

  }

  #axis labeling, be careful b/c of coord_flip
  if(is.null(ylab)) ylab = ""
  if(is.null(xlab)) xlab = "Posterior Estimate"

  final = ggplot(plot.df, aes(x=names, y=c.point)) +
    geom_errorbar(ymin=c.lower, ymax=c.upper, width=0, col="grey60") +
    geom_point() +
    geom_point(data = plot.df, aes(x=names, y=c.upper), size=-1) + #for correct scaling
    geom_point(data = plot.df, aes(x=names, y=c.lower), size=-1) + #for correct scaling
    theme_bw()   +
    geom_hline(yintercept = 0,lty=2) +
    coord_flip() +
    xlab(ylab)   +
    ylab(xlab)

  return(final)

}



#' @name plot.UPG.MNL
#' @title Coefficient plots for UPG.MNL objects
#'
#' @description \code{plot} generates plots from \code{UPG.MNL} objects using \code{ggplot2}. Coefficient plots show point estimates for all coefficients in all groups except the baseline as well as their credible intervals.
#'
#' @param x an object of class \code{UPG.MNL}.
#' @param sort a logical variable indicating whether the plotted coefficients should be sorted according to average effect sizes across groups. Default is FALSE.
#' @param names a character vector indicating names for the variables used in the plots.
#' @param groups  a character vector indicating names for the groups excluding the baseline. The group names must correspond to the ordering in the dependent variable used for estimation.
#' @param xlab  a character vector of length 1 indicating a title for the x-axis.
#' @param ylab  a character vector of length 1 indicating a title for the y-axis.
#' @param q a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#' @param include can be used to plot only a subset of variables. Specify the columns of X that should be kept in the plot. See examples for further information.
#' @param ... other plot parameters.
#'
#' @return Returns a ggplot2 object.
#'
#' @seealso
#' \code{\link{summary.UPG.MNL}} to summarize a \code{UPG.MNL} object and create tables.
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
#' # plot the results and sort coefficients by average effect size
#' plot(results.mnl, sort = TRUE)
#'
#' # plot only variables 1 and 3 with custom group and variable names
#' # also, customize credible intervals and axis labels
#' plot(results.mnl,
#'      include  = c(1,3),
#'      names    = c("Custom 1", "Custom 2"),
#'      groups   = c("Alpha", "Beta"),
#'      q        = c(0.1, 0.9),
#'      xlab     = c("Custom X"),
#'      ylab     = c("Custom Y"))
#'}
#' @method  plot UPG.MNL
#'
#'@export
plot.UPG.MNL    = function(x         = NULL,
                           ...,
                           sort      = FALSE,           # sort coefficients by average effect size
                           names     = NULL,            # provide names for variables, alternatively
                           groups    = NULL,            # provide names for groups except baseline
                           xlab      = NULL,            # provide x axis label
                           ylab      = NULL,            # provide y axis label
                           q         = c(0.025, 0.975), # credible intervals
                           include   = NULL             # which variables to include? default:all (numeric vector)
){



  if(is.null(include)) include = 1:ncol(x$inputs$X)
  if(is.null(names))   names = colnames(x$inputs$X[,include,drop=F])
  if(is.null(names))   names = paste0("Variable", 1:ncol(x$inputs$X[,include,drop=F]))
  if(length(names) != length(include)) stop("Number of provided variable names does not match number of included variables.")


  #create some global variables such that r cmd check is happy
  variable = iter = value = NULL

  c.point  = apply(x$posterior$beta[,include,,drop=F], c(2,3), mean)
  c.upper  = apply(x$posterior$beta[,include,,drop=F], c(2,3), quantile, q[2])
  c.lower  = apply(x$posterior$beta[,include,,drop=F], c(2,3), quantile, q[1])


  if(nrow(c.point) == 1){c.upper = matrix(c.upper,nrow=1); c.lower = matrix(c.lower,nrow=1)}

  if(is.null(groups)) groups = paste(x$posterior$groups[-ncol(c.point)])
  if(length(groups) != ncol(c.point)-1) stop("Wrong number of group names supplied. Need K-1 names where K is the number of choices.")


  #kick baseline
  c.upper = c.upper[,-ncol(c.upper),drop=F]
  c.lower = c.lower[,-ncol(c.lower),drop=F]
  c.point = c.point[,-ncol(c.point),drop=F]

  #add variable names
  c.upper = data.frame(c.upper, names)
  c.lower = data.frame(c.lower, names)
  c.point = data.frame(c.point, names)

  #add group names
  colnames(c.upper) = colnames(c.lower) = colnames(c.point) = c(groups, "names")

  #add measurement
  c.upper$measure   = "c.upper"
  c.lower$measure   = "c.lower"
  c.point$measure   = "c.point"

  plot.df = rbind(c.upper,c.lower,c.point)
  plot.df = reshape2::melt(plot.df, id.vars = c("names","measure"))
  plot.df = reshape2::dcast(plot.df, "names + variable ~ measure")

  #sorting (bit more complicated in MNL, I use average point estimate over all groups)
  if(length(names)>1){
    if(sort){

      average       = aggregate(plot.df$c.point, by=list(plot.df$names), FUN=mean)
      lvls          = unique(plot.df$names)
      plot.df$names = factor(plot.df$names, levels = lvls[order(average$x)])


    }
  }

  if(!sort){

    # plot in order of appearance in X
    plot.df$names = factor(plot.df$names, levels = rev(names))

  }

  #axis labeling, take care b/c of coord_flip
  if(is.null(ylab)) ylab = ""
  if(is.null(xlab)) xlab = "Posterior Estimate"


   final =  ggplot(plot.df, aes(x=names, y=c.point, shape = variable, group = variable)) +
             geom_errorbar(aes(x=names, ymin=c.lower, ymax=c.upper, group = variable),
                               position = position_dodge(width=.5),
                               width=0, col="grey60") +
             geom_point(position = position_dodge(width=.5)) +
             theme_bw()   +
             geom_hline(yintercept = 0,lty=2) +
             coord_flip() +
             xlab(ylab)     +
             ylab(xlab) +
             theme(legend.position = "bottom",
                   legend.title = element_blank())

   return(final)
}


#' @name plot.UPG.Binomial
#' @title Coefficient plots for UPG.Binomial objects
#'
#' @description \code{plot} generates plots from \code{UPG.Binomial} objects using \code{ggplot2}. Coefficient plots show point estimates for all coefficients as well as their credible intervals.
#'
#' @param x an object of class \code{UPG.Binomial}.
#' @param sort a logical variable indicating whether the plotted coefficients should be sorted according to effect sizes. Default is FALSE.
#' @param names a character vector indicating names for the variables used in the plots.
#' @param xlab  a character vector of length 1 indicating a title for the x-axis.
#' @param ylab  a character vector of length 1 indicating a title for the y-axis.
#' @param q a numerical vector of length two providing the posterior quantiles to be extracted. Default are 0.025 and 0.975 quantiles.
#' @param include can be used to plot only a subset of variables. Specify the columns of X that should be kept in the plot. See examples for further information.
#' @param ... other plot parameters.
#'
#' @return Returns a ggplot2 object.
#'
#' @seealso
#' \code{\link{summary.UPG.Binomial}} to summarize a \code{UPG.Binomial} object and create tables.
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
#' # plot the results and sort coefficients by effect size
#' plot(results.binomial, sort = TRUE)
#'
#' # plot only variables 1 and 3 with custom names, credible intervals and axis labels
#' plot(results.binomial,
#'      include  = c(1,3),
#'      names    = c("Custom 1", "Custom 2"),
#'      q        = c(0.1, 0.9),
#'      xlab     = c("Custom X"),
#'      ylab     = c("Custom Y"))
#'}
#' @method  plot UPG.Binomial
#'
#'@export
plot.UPG.Binomial  = function(x         = NULL,
                              ...,
                              sort      = FALSE,           # sort coefficients by average effect size
                              names     = NULL,            # provide names for variables, alternatively
                              xlab      = NULL,            # provide x axis label
                              ylab      = NULL,            # provide y axis label
                              q         = c(0.025, 0.975), # credible intervals
                              include   = NULL             # which variables to include? default:all (numeric vector)
){



  if(is.null(include)) include = 1:ncol(x$inputs$X)
  if(is.null(names))   names = colnames(x$inputs$X[,include,drop=F])
  if(is.null(names))   names = paste0("Variable", 1:ncol(x$inputs$X[,include,drop=F]))
  if(length(names) != length(include)) stop("Number of provided variable names does not match number of included variables.")

  #create some global variables such that r cmd check is happy
  variable = iter = value = NULL

  c.point  = apply(x$posterior$beta[,include,drop=F], 2, mean)
  c.upper  = apply(x$posterior$beta[,include,drop=F], 2, quantile, q[2])
  c.lower  = apply(x$posterior$beta[,include,drop=F], 2, quantile, q[1])

  plot.df  = data.frame(c.point,c.upper,c.lower,names)

  if(nrow(plot.df) > 1 & sort){

    plot.df$names = factor(plot.df$names, levels = plot.df$names[order(plot.df$c.point)])

  }

  if(!sort){

    # plot in order of appearance in X
    plot.df$names = factor(plot.df$names, levels = rev(plot.df$names))

  }

  #axis labeling, be careful b/c of coord_flip
  if(is.null(ylab)) ylab = ""
  if(is.null(xlab)) xlab = "Posterior Estimate"

  final = ggplot(plot.df, aes(x=names, y=c.point)) +
    geom_errorbar(ymin=c.lower, ymax=c.upper, width=0, col="grey60") +
    geom_point() +
    geom_point(data = plot.df, aes(x=names, y=c.upper), size=-1) + #for correct scaling
    geom_point(data = plot.df, aes(x=names, y=c.lower), size=-1) + #for correct scaling
    theme_bw()   +
    geom_hline(yintercept = 0,lty=2) +
    coord_flip() +
    xlab(ylab)   +
    ylab(xlab)

  return(final)

}
