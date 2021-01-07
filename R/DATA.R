#' Students program choices.
#'
#' A dataset containing the choice among general program, vocational program and academic program for 200 high school students as well as some explanatory
#' variables.
#'
#' @format A data frame with 200 rows and 5 variables:
#' \describe{
#'   \item{program}{a vector of program choices}
#'   \item{intercept}{an intercept}
#'   \item{female}{dummy for female students}
#'   \item{ses}{socioeconomic status, 1 is lowest}
#'   \item{write}{writing score of student}
#' }
#' @source Original dataset is known as 'hsbdemo' and has been sourced from \url{https://stats.idre.ucla.edu/stat/data/hsbdemo.dta}.
"program"

#' Grouped Titanic survival data.
#'
#' A dataset containing the number of survivals and the total number of persons by
#' passenger class, age group and gender.
#'
#' @format A data frame with 78 rows and 6 variables:
#' \describe{
#'   \item{survived}{number of passengers that survived}
#'   \item{total}{number of total passengers}
#'   \item{intercept}{an intercept}
#'   \item{pclass}{passenger class (1 highest, 3 lowest)}
#'   \item{female}{dummy for females}
#'   \item{age.group}{age group indicator (0-5yrs, 5-10yrs, ...)}
#' }
#' @source Data originally sourced from \url{https://web.stanford.edu/class/archive/cs/cs109/cs109.1166/stuff/titanic.csv}. See also
#' \url{https://towardsdatascience.com/the-binomial-regression-model-everything-you-need-to-know-5216f1a483d3}.
"titanic"

#' Female labor force participation data.
#'
#' A dataset containing socio-economic characteristics as well as a labor
#' force participation dummy for 753 married women from the panel study
#' of income dynamics.
#'
#' @format A data frame with 753 rows and 9 variables:
#' \describe{
#'   \item{lfp}{binary indicator for participating in the labor force (=1) or not (=0)}
#'   \item{intercept}{intercept}
#'   \item{k5}{number of children 5 years old or younger}
#'   \item{k618}{number of children 6 to 18 years old}
#'   \item{age}{in years}
#'   \item{wc}{wife went to college dummy}
#'   \item{hc}{husband went to college dummy}
#'   \item{lwg}{log expected wage rate; for women in the labor force, the actual wage rate; for women not in the labor force, an imputed value based on the regression of \code{lwg} on the other variables}
#'   \item{inc}{family income exclusive of wife's income}
#' }
#' @source Data taken from 'carData' package. Also known as the 'Mroz' dataset. Mroz, T. A. (1987) The sensitivity of an empirical model of married women's hours of work to economic and statistical assumptions. Econometrica 55, 765-799.
"lfp"
