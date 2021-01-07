#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace Rcpp;

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
double aterm_l(int n, double x, double t)
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
double randinvg_l(double mu)
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
double exprnd_l(double mu)
{
  return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

// Sample truncated gamma random variates
// Ref: Chung, Y.: Simulation of truncated gamma variables
// Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
double truncgamma_l()
{
  double c = MATH_PI_2;
  double X, gX;

  bool done = false;
  while(!done)
  {
    X = exprnd_l(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);

    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }

  return X;
}

// Sample truncated inverse Gaussian random variates
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
double tinvgauss_l(double z, double t)
{
  double X, u;
  double mu = 1.0/z;

  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / truncgamma_l();

      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = randinvg_l(mu);
    }
  }
  return X;
}






// Sample PG(1,z)
// Based on Algorithm 6 in PhD thesis of Jesse Bennett Windle, 2013
// URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
double samplepg_l(double z)
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
      X = t + exprnd_l(1.0)/K;
    }
    else {
      // truncated Inverse Gaussian
      X = tinvgauss_l(z, t);
    }

    // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    int i = 1;
    double Sn = aterm_l(0, X, t);
    double U = R::runif(0.0,1.0) * Sn;
    int asgn = -1;
    bool even = false;

    while(1)
    {
      Sn = Sn + asgn * aterm_l(i, X, t);

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





NumericVector rcpp_pgdraw_l(NumericVector b, NumericVector c)
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
      y[i] += samplepg_l(c[i]);
    }
  }

  return y;
}


arma::vec Fel_(arma::vec x){

  NumericVector x1 = as<NumericVector>(wrap(x));
  arma::vec logged = log(x1)-log(1-x1);
  return logged;

  }


arma::vec Fel(arma::vec x){

  NumericVector x1 = as<NumericVector>(wrap(x));
  arma::vec exped = exp(x1)/(1+exp(x1));
  return exped;

  }


//' @name sample_logit
//' @noRd
// [[Rcpp::export]]
List sample_logit( arma::mat y_matrix,
                   arma::mat y_categ,
                   arma::mat X,
                   arma::mat beta_old,
                   int       nsave,
                   int       nburn,
                   double    A0,
                   double    d0,
                   double    D0,
                   double    G0,
                   int       BOOST,
                   int       verbose){

  // GIBBS SAMPLER SETUP

  int N    = y_matrix.n_rows;
  int P    = X.n_cols;
  // int K    = y_matrix.n_cols;
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

  arma::cube beta_store(nsave, P, 1); beta_store.fill(0);
  arma::mat  gamma_store(nsave, 1);   gamma_store.fill(0);
  arma::mat  delta_store(nsave, 1);   delta_store.fill(0);
  arma::cube y_lat_store(nsave,N,1);  y_lat_store.fill(0);

  // STARTING VALUES

  arma::mat beta_draw = beta_old;

  // PROGRESS BAR
  Progress p(ntot, verbose==1);

  // GIBBS SAMPLER

  for(int irep = 0; irep < ntot; irep++){

    // SAMPLE LATENT UTILITIES

    arma::mat U_i    = Rcpp::runif(N, 0, 1);
    arma::mat lambda = X * beta_draw;
    arma::mat y_d    = lambda + Fel_(y_matrix.col(0) + U_i % (1 - y_matrix.col(0) - (1 - Fel(-lambda))));


    // SAMPLE SCALING FACTORS

      NumericVector b(N); b.fill(2);
      arma::vec     c  = y_d - lambda;
      NumericVector c1 = as<NumericVector>(wrap(c));
      NumericVector c2 = Rcpp::abs(c1);
      arma::vec omega  = 1 / rcpp_pgdraw_l(b,c2);




    // DRAW DELTA.STAR

    double delta_star = 1/R::rgamma(d0, 1/D0);

    // DRAW GAMMA.STAR

    double g1 = - beta_draw.at(0,0) * (G0/(G0+A0));
    double gamma_star = R::rnorm(sqrt(delta_star) * g1, sqrt(delta_star * G1));


    if(BOOST==0){
      gamma_star = 0;
      delta_star = 1;
      }

    // SWITCH TO UNIDENTIFIED MODEL
    arma::mat y_d_tilde = y_d * sqrt(delta_star) + gamma_star;

    // DRAW GAMMA
    double u = R::runif(0,1);

    arma::vec y1  = y_matrix.col(0);
    arma::vec yd1 = y_d_tilde.rows(find(y1 == 1));
    arma::vec yd0 = y_d_tilde.rows(find(y1 == 0));

    double O = 0;
    double L = 0;

    if(sum(y1) == 0){O = R_PosInf;};
    if(sum(y1) != 0){O = min(yd1);};
    if(sum(y1) == N){L = R_NegInf;};
    if(sum(y1) != N){L = max(yd0);};


    double O_ = O/sqrt(G0star);
    double L_ = L/sqrt(G0star);

    double prob  = u  * R::pt(O_, 2*d0, 1,0) + (1 - u) * R::pt(L_, 2*d0, 1, 0);
    double gamma = sqrt(G0star) * R::qt(prob, 2*d0, 1, 0);

    // DRAW COEFFICIENTS

       // TRANSFORM TO 'UNIT VARIANCE MODEL'

       arma::mat y_star = y_d_tilde / sqrt(omega);
       arma::mat X_star = X;

       for(int pp = 0; pp < P; pp++){

         X_star.col(pp) = X_star.col(pp) / sqrt(omega);

       }

       // DRAW Gn AND gn IN UNIDENTIFIED MODEL

       arma::mat G_n_inv = A0_inv + X_star.t() * X_star;
       arma::mat G_n     = inv(G_n_inv);
       arma::mat g_n     = G_n * (X_star.t()   * y_star);

       // DRAW DELTA
       arma::vec resids = y_d_tilde - X * g_n;
       double SSR       = sum(pow(resids,2) / omega);
       arma::vec Dn0    = D0 + SSR/2 + pow(gamma,2) / (2*G0) + (g_n.t() * A0_inv * g_n)/2;
       double Dn        = Dn0(0);

       double delta = 1/R::rgamma(dn, 1/Dn);


     if(BOOST == 0){
         gamma = 0;
         delta = 1;
      }


       // DRAW FROM MVN

       arma::mat      sig =  delta * G_n;
       arma::vec      mu  =  g_n;

       arma::vec rand_n   = rnorm(sig.n_cols);
       arma::mat sig_cho  = arma::chol(sig,"lower");
       arma::vec draw     = mu + sig_cho * rand_n;

       beta_draw = draw;

   // TRANSFORM BACK
   arma::mat gamma_trans(P,1); gamma_trans.fill(0);
   gamma_trans.at(0,0) = gamma;

   beta_draw = beta_draw - gamma_trans;
   beta_draw = beta_draw / sqrt(delta);

    // STORE DRAWS

       if(irep >= nburn){

         beta_store.row(irep  - nburn)  = beta_draw;
         y_lat_store.row(irep - nburn)  = y_d;
         gamma_store.row(irep - nburn)  = gamma;
         delta_store.row(irep - nburn)  = delta;

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




















