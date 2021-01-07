#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>


// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;





/*
* This file implements the Polya-gamma sampler PG(1,z).
* This is a C++ implementation of Algorithm 6 in PhD thesis of Jesse
* Bennett Windle, 2013
* URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
*
* References:
*
*   Jesse Bennett Windle
*   Forecasting High-Dimensional, Time-Varying Variance-Covariance Matrices
*   with High-Frequency Data and Sampling Polya-Gamma Random Variates for
*   Posterior Distributions Derived from Logistic Likelihoods
*   PhD Thesis, 2013
*
*   Damien, P. & Walker, S. G. Sampling Truncated Normal, Beta, and Gamma Densities
*   Journal of Computational and Graphical Statistics, 2001, 10, 206-215
*
*   Chung, Y.: Simulation of truncated gamma variables
*   Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
*
* (c) Copyright Enes Makalic and Daniel F Schmidt, 2018
*/

// Mathematical constants computed using Wolfram Alpha
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392


// Function a_n(x) defined in equations (12) and (13) of
// Bayesian inference for logistic models using Polya-Gamma latent variables
// Nicholas G. Polson, James G. Scott, Jesse Windle
// arXiv:1205.0310
//
// Also found in the PhD thesis of Windle (2013) in equations
// (2.14) and (2.15), page 24
double aterm_m(int n, double x, double t)
{
  double f = 0;
  if(x <= t) {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  }
  else {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }
  return (double)exp(f);
}

// Generate inverse gaussian random variates
double randinvg_m(double mu)
{
  // sampling
  double u = R::rnorm(0.0,1.0);
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );

  if(R::runif(0.0,1.0) > mu /(mu+out)) {
    out = mu*mu / out;
  }
  return out;
}

// Generate exponential distribution random variates
double exprnd_m(double mu)
{
  return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

// Sample truncated gamma random variates
// Ref: Chung, Y.: Simulation of truncated gamma variables
// Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
double truncgamma_m()
{
  double c = MATH_PI_2;
  double X, gX;

  bool done = false;
  while(!done)
  {
    X = exprnd_m(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);

    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }

  return X;
}

// Sample truncated inverse Gaussian random variates
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
double tinvgauss_m(double z, double t)
{
  double X, u;
  double mu = 1.0/z;

  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / truncgamma_m();

      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = randinvg_m(mu);
    }
  }
  return X;
}






// Sample PG(1,z)
// Based on Algorithm 6 in PhD thesis of Jesse Bennett Windle, 2013
// URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
double samplepg_m(double z)
{
  //  PG(b, z) = 0.25 * J*(b, z/2)
  z = (double)std::fabs((double)z) * 0.5;

  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
  double t = MATH_2_PI;

  // Compute p, q and the ratio q / (q + p)
  // (derived from scratch; derivation is not in the original paper)
  double K = z*z/2.0 + MATH_PI2/8.0;
  double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
  double logK = (double)std::log(K);
  double Kt = K * t;
  double w = (double)std::sqrt(MATH_PI_2);

  double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
  double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
  double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
  double ratio = 1.0 / (1.0 + p_over_q);

  double u, X;

  // Main sampling loop; page 130 of the Windle PhD thesis
  while(1)
  {
    // Step 1: Sample X ? g(x|z)
    u = R::runif(0.0,1.0);
    if(u < ratio) {
      // truncated exponential
      X = t + exprnd_m(1.0)/K;
    }
    else {
      // truncated Inverse Gaussian
      X = tinvgauss_m(z, t);
    }

    // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    int i = 1;
    double Sn = aterm_m(0, X, t);
    double U = R::runif(0.0,1.0) * Sn;
    int asgn = -1;
    bool even = false;

    while(1)
    {
      Sn = Sn + asgn * aterm_m(i, X, t);

      // Accept if n is odd
      if(!even && (U <= Sn)) {
        X = X * 0.25;
        return X;
      }

      // Return to step 1 if n is even
      if(even && (U > Sn)) {
        break;
      }

      even = !even;
      asgn = -asgn;
      i++;
    }
  }
  return X;
}





NumericVector rcpp_pgdraw_m(NumericVector b, NumericVector c)
{
  int m = b.size();
  int n = c.size();
  NumericVector y(n);

  // Setup
  int i, j, bi = 1;
  if (m == 1)
  {
    bi = b[0];
  }

  // Sample
  for (i = 0; i < n; i++)
  {
    if (m > 1)
    {
      bi = b[i];
    }

    // Sample
    y[i] = 0;
    for (j = 0; j < (int)bi; j++)
    {
      y[i] += samplepg_m(c[i]);
    }
  }

  return y;
}


arma::mat RcppRowmaxs(arma::mat x){

  int N = x.n_rows;
  arma::mat rm(N,1);

  for(int nn = 0; nn < N; nn++){

    rm(nn) = max(x.row(nn));

  }

  return(rm);

}

arma::mat sample_utils(arma::mat lambda_star,
                       arma::mat lambda_k,
                       arma::mat lambda_a,
                       arma::mat y_matrix,
                       int categ,
                       arma::mat y_categ,
                       int K0){

  int N = lambda_star.n_rows;

  arma::mat ut_matrix(N,3); ut_matrix.fill(0);
  arma::mat Ui     = Rcpp::runif(N, 0, 1);
  arma::mat Vki    = Rcpp::runif(N, 0, 1);
  arma::mat V0i    = Rcpp::runif(N, 0, 1);
  arma::mat Vai    = Rcpp::runif(N, 0, 1);

  // INDICATORS
  arma::vec ind_categ(N,1); ind_categ.fill(0);
  arma::vec ind_other(N,1); ind_other.fill(0);
  arma::vec ind_basel(N,1); ind_basel.fill(0);

  LogicalVector lc  = as<LogicalVector>(wrap( y_categ  != categ));
  LogicalVector lo1 = as<LogicalVector>(wrap((y_categ  == categ)));
  LogicalVector lo2 = as<LogicalVector>(wrap((y_categ  == (K0-1))));
  LogicalVector lo  = lo1 | lo2;
  LogicalVector lb  = as<LogicalVector>(wrap( y_categ  != (K0-1)));

  for(int nn = 0; nn < N; nn++){
    if(lc(nn)){ind_categ(nn) = 1;}
    if(lo(nn)){ind_other(nn) = 1;}
    if(lb(nn)){ind_basel(nn) = 1;}
  }

  // FIRST PART OF DIFFERENCES
  arma::mat utilities = - log(Ui) % (1/lambda_star);

  // INDIVIDUAL IS CATEG
  ut_matrix.col(0) = utilities - (log(Vki) % (1/lambda_k)) % ind_categ;

  // INDIVIDUAL IS OTHER
  ut_matrix.col(1) = utilities - (log(Vai) % (1/lambda_a)) % ind_other;

  // INDIVIDUAL IS BASELINE
  ut_matrix.col(2) = utilities - log(V0i) % ind_basel;

  arma::mat utilities_final =  - log(ut_matrix);


  return(utilities_final);

}


//' @name sample_mnl
//' @noRd
// [[Rcpp::export]]
List sample_mnl(   arma::mat  y_matrix,
                   arma::mat  y_categ,
                   arma::mat  X,
                   arma::mat  beta_old,
                   int        nsave,
                   int        nburn,
                   double     A0,
                   double     d0,
                   double     D0,
                   double     G0,
                   int        BOOST,
                   int        verbose){



  // GIBBS SAMPLER SETUP
  int N    = y_matrix.n_rows;
  int P    = X.n_cols;
  int K    = y_matrix.n_cols;
  int ntot = nburn + nsave;

  // PRIOR SETUP
  arma::vec  priorvecs(P); priorvecs.fill(A0);

  if(BOOST == 0){
    priorvecs(0) = G0;
  }

  arma::mat  A0p    = diagmat(priorvecs);
  arma::mat  A0_inv = inv(A0p);

  double G0star  = D0  * G0/d0;
  double G1      = (G0 * A0) / (G0 + A0);
  double dn      = d0 + ((N+1)/2);

  // STORAGE
  arma::cube beta_store(nsave,  P, K);  beta_store.fill(0);
  arma::mat  gamma_store(nsave, K);     gamma_store.fill(0);
  arma::mat  delta_store(nsave, K);     delta_store.fill(0);
  arma::cube y_lat_store(nsave, N, 3);  y_lat_store.fill(0);

  // STARTING VALUES
  arma::mat y_d(N,3);
  arma::mat beta_draw = beta_old;
  arma::vec delta(K); delta.fill(1);
  arma::vec gamma(K); gamma.fill(0);

  // PROGRESS BAR
  Progress p(ntot, verbose==1);

  arma::vec zeroes(N); zeroes.fill(0);

  // GIBBS SAMPLER

  for(int irep = 0; irep < ntot; irep++){

    for(int kk = 0; kk < (K-1); kk++){

    // SAMPLE LATENT UTILITIES

      arma::mat lambda = exp(X * beta_draw);
      arma::mat lambda_star(N,1);

      lambda_star = sum(lambda,1);

      arma::mat  lambda_k    = lambda.col(kk);
      arma::mat  lambda_d    = lambda;
      lambda_d.col(kk)       = zeroes;
      lambda_d.col(K-1)      = zeroes;
      arma::mat  lambda_a    = sum(lambda_d, 1);

      y_d         = sample_utils(lambda_star,
                                 lambda_k,
                                 lambda_a,
                                 y_matrix,
                                 kk,
                                 y_categ,
                                 K);


    // SAMPLE SCALING FACTORS

      NumericVector b(N); b.fill(2);
      arma::vec     c  = y_d.col(0) - y_d.col(2) - log(lambda_k);
      NumericVector c1 = as<NumericVector>(wrap(c));
      NumericVector c2 = Rcpp::abs(c1);
      arma::vec omega  = 1 / rcpp_pgdraw_m(b,c2);

    // DRAW DELTA.STAR

    double delta_star = 1/R::rgamma(d0, 1/D0);

    // DRAW GAMMA.STAR

    double g1 = - beta_draw.at(0,kk) * (G0/(G0+A0));
    double gamma_star = R::rnorm(sqrt(delta_star) * g1, sqrt(delta_star * G1));


    if(BOOST==0){
      gamma_star = 0;
      delta_star = 1;
      }

    // SWITCH TO UNIDENTIFIED MODEL
    arma::mat y_d_tilde(N,3);
    y_d_tilde.col(0) = y_d.col(0) * sqrt(delta_star) + gamma_star;
    y_d_tilde.col(1) = y_d.col(1) * sqrt(delta_star);
    y_d_tilde.col(2) = y_d.col(2) * sqrt(delta_star);

    // DRAW GAMMA
    double u = R::runif(0,1);

    // DO A LOT OF THINGS FOR THE IF-ELSE STATEMENTS
    arma::vec yk      = y_matrix.col(kk);
    arma::vec yb      = y_matrix.col(K-1);

    arma::vec y1y3    = y_d_tilde.col(0) - y_d_tilde.col(2);
    arma::vec y1y31   = y1y3(find(yb==1));

    arma::vec y1y2    = y_d_tilde.col(0) - y_d_tilde.col(1);

    LogicalVector yc1  = as<LogicalVector>(wrap((yb  == 0)));
    LogicalVector yc2  = as<LogicalVector>(wrap((yk  == 0)));
    LogicalVector ycc  = yc1 & yc2;

    arma::vec     ycv(N); ycv.fill(0);

    for(int nn = 0; nn < N; nn++){
      if(ycc(nn)){ycv(nn) = 1;}
    }

    arma::vec y1y20  = y1y2(find(ycv == 1));

    arma::mat y1          = y_d_tilde.col(0);
    arma::mat y23(N,2); y23.col(0) = y_d_tilde.col(1); y23.col(1) = y_d_tilde.col(2);
    arma::mat rmy23       = RcppRowmaxs(y23);
    arma::vec y1rmy23     = y1 - rmy23;
    arma::vec y1rmy231    = y1rmy23(find(yk==1));

    double max1 = 0;
    double max2 = 0;
    double O = 0;
    double L = 0;

    if(sum(yk)  == 0)    {O    = R_PosInf;};
    if(sum(yk)  != 0)    {O    = min(y1rmy231);};

    if(sum(yb == 1) == 0){max1 = R_NegInf;};
    if(sum(yb == 1) != 0){max1 = max(y1y31);};


    LogicalVector l1 = as<LogicalVector>(wrap(yb  == 0));
    LogicalVector l2 = as<LogicalVector>(wrap(yk  == 0));
    LogicalVector l  = l1 & l2;

    if(sum(l) == 0){max2 = R_NegInf;};
    if(sum(l) != 0){max2 = max(y1y20);};

    arma::vec choice(2); choice(0) = max1; choice(1) = max2;
    L = max(choice);

    double O_ = O/sqrt(G0star);
    double L_ = L/sqrt(G0star);

    double prob  = u  * R::pt(O_, 2*d0, 1,0) + (1 - u) * R::pt(L_, 2*d0, 1, 0);
    gamma(kk)    = sqrt(G0star) * R::qt(prob, 2*d0, 1, 0);

    // DRAW COEFFICIENTS

       // TRANSFORM TO 'UNIT VARIANCE MODEL'

       arma::mat y_dep  = y_d_tilde.col(0) - y_d_tilde.col(2);
       arma::mat y_star = y_dep / sqrt(omega);
       arma::mat X_star = X;

       for(int pp = 0; pp < P; pp++){

         X_star.col(pp) = X_star.col(pp) / sqrt(omega);

       }

       // DRAW Gn AND gn IN UNIDENTIFIED MODEL
       arma::mat G_n_inv = A0_inv + X_star.t() * X_star;
       arma::mat G_n     = inv(G_n_inv);
       arma::mat g_n     = G_n * (X_star.t()   * y_star);

       // DRAW DELTA
       arma::vec resids = y_dep - X * g_n;
       double SSR       = sum(pow(resids,2) / omega);
       arma::vec Dn0    = D0 + SSR/2 + pow(gamma(kk),2) / (2*G0) + (g_n.t() * A0_inv * g_n)/2;
       double Dn        = Dn0(0);

       delta(kk) = 1/R::rgamma(dn, 1/Dn);


     if(BOOST == 0){
         gamma(kk) = 0;
         delta(kk) = 1;
      }


       // DRAW FROM MVN

       arma::mat      sig =  delta(kk) * G_n;
       arma::vec      mu  =  g_n;

       arma::vec rand_n   = rnorm(sig.n_cols);
       arma::mat sig_cho  = arma::chol(sig,"lower");
       arma::vec draw     = mu + sig_cho * rand_n;

       beta_draw.col(kk)  = draw;

   // TRANSFORM BACK
   arma::mat gamma_trans(P,1); gamma_trans.fill(0);
   gamma_trans.at(0,0) = gamma(kk);

   beta_draw.col(kk) = beta_draw.col(kk) - gamma_trans;
   beta_draw.col(kk) = beta_draw.col(kk) / sqrt(delta(kk));

    }

    // STORE DRAWS

       if(irep >= nburn){

         beta_store.row(irep  - nburn)  = beta_draw;
         y_lat_store.row(irep - nburn)  = y_d;
         gamma_store.row(irep - nburn)  = gamma.t();
         delta_store.row(irep - nburn)  = delta.t();

        }

    // PROGRESS BAR

    if(verbose == 1) p.increment(); // update progress

    // CHECK USER INTERRUPTION

    if(irep % 20 == 0) Rcpp::checkUserInterrupt();


}





return List::create(Named("beta")     = beta_store,
                    Named("y.latent") = y_lat_store,
                    Named("gamma")    = gamma_store,
                    Named("delta")    = delta_store);


}




















