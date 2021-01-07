// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <progress.hpp>
#include <progress_bar.hpp>

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
double aterm_b(int n, double x, double t)
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
double randinvg_b(double mu)
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
double exprnd_b(double mu)
{
  return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

// Sample truncated gamma random variates
// Ref: Chung, Y.: Simulation of truncated gamma variables
// Korean Journal of Computational & Applied Mathematics, 1998, 5, 601-610
double truncgamma_b()
{
  double c = MATH_PI_2;
  double X, gX;

  bool done = false;
  while(!done)
  {
    X = exprnd_b(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);

    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }

  return X;
}

// Sample truncated inverse Gaussian random variates
// Algorithm 4 in the Windle (2013) PhD thesis, page 129
double tinvgauss_b(double z, double t)
{
  double X, u;
  double mu = 1.0/z;

  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / truncgamma_b();

      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = randinvg_b(mu);
    }
  }
  return X;
}






// Sample PG(1,z)
// Based on Algorithm 6 in PhD thesis of Jesse Bennett Windle, 2013
// URL: https://repositories.lib.utexas.edu/bitstream/handle/2152/21842/WINDLE-DISSERTATION-2013.pdf?sequence=1
double samplepg_b(double z)
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
      X = t + exprnd_b(1.0)/K;
    }
    else {
      // truncated Inverse Gaussian
      X = tinvgauss_b(z, t);
    }

    // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    int i = 1;
    double Sn = aterm_b(0, X, t);
    double U = R::runif(0.0,1.0) * Sn;
    int asgn = -1;
    bool even = false;

    while(1)
    {
      Sn = Sn + asgn * aterm_b(i, X, t);

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





NumericVector rcpp_pgdraw_b(NumericVector b, NumericVector c)
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
      y[i] += samplepg_b(c[i]);
    }
  }

  return y;
}

double d_sq(double bi, double BI){
  return(pow(BI,2) / (4*pow((bi+1),2)));
}

double d_ig(double bi, double dI, double DI){
  return(DI/(dI-bi));
}

double d_m(double DI, double dI, double BI){
  return(16*pow(DI,2) / pow((BI + sqrt(pow(BI,2) + 16 * DI * (dI+1))),2));
}

double I_p(double DI, double dI, double BI, double dm){
  return(- (sqrt(pow(BI,2) + 16*DI*(dI+1))) / (4*pow(dm,5/2)));
}


IntegerVector seq_cpp(int lo, int hi) {
  int n = hi - lo + 1;

  // Create a new integer vector, sequence, of size n
  IntegerVector sequence(n);

  for(int i = 0; i < n; i++) {
    // Set the ith element of sequence to lo plus i
    sequence(i) = lo + i;
  }

  return sequence;
}

NumericVector vecpow(const NumericVector base, const NumericVector exp) {

  NumericVector out(base.size());

  std::transform(base.begin(), base.end(),
                 exp.begin(), out.begin(), ::powf);
  return out;
}

mat submat_arma(arma::mat X, arma::colvec T) {
  mat y = X.rows(find(T == 1));
  return y;
}


//' @name sample_binom
//' @noRd
// [[Rcpp::export]]

List sample_binom( arma::mat y_matrix,
                   arma::mat Ni,
                   arma::mat X,
                   arma::mat beta_old,
                   int       nsave,
                   int       nburn,
                   double    A0,
                   double    d0,
                   double    D0,
                   double    G0,
                   int       BOOST,
                   int       verbose,
                   int       prior_ig){

  // GIBBS SAMPLER SETUP
  int N    = y_matrix.n_rows;
  int P    = X.n_cols;
  int ntot = nburn + nsave;

  // INDICATORS FOR BINOMIAL MODEL
  arma::uvec     Wix  = find(y_matrix > 0);
  arma::uvec     Vix  = find(y_matrix < Ni);
  LogicalVector Windl = as<LogicalVector>(wrap(y_matrix > 0));
  LogicalVector Vindl = as<LogicalVector>(wrap(y_matrix < Ni));
  arma::vec Wind(N);
  arma::vec Vind(N);

  for(int nn = 0; nn < N; nn++){
    if(Windl(nn)){Wind(nn) = 1;}else{Wind(nn) = 0;}
    if(Vindl(nn)){Vind(nn) = 1;}else{Vind(nn) = 0;}
  }

  int NW = sum(Wind);
  int NV = sum(Vind);


  // KAPPA
  arma::mat kappa_w = (1  - y_matrix) / 2;
  arma::mat kappa_v = (Ni - y_matrix - 1) / 2;
  arma::mat kappa_ww = submat_arma(kappa_w, Wind);
  arma::mat kappa_vv = submat_arma(kappa_v, Vind);

  // PRIOR SETUP

  arma::vec  priorvecs(P); priorvecs.fill(A0);

  if(BOOST == 0){
    priorvecs(0) = G0;
  }

  arma::mat  A0p    = diagmat(priorvecs);
  arma::mat  A0_inv = inv(A0p);

  double G0star  = D0  * G0/d0;
  double G1      = (G0 * A0) / (G0 + A0);

  // STORAGE

  arma::cube  beta_store(nsave, P, 1); beta_store.fill(0);
  arma::mat   gamma_store(nsave, 1);   gamma_store.fill(0);
  arma::mat   delta_store(nsave, 1);   delta_store.fill(0);
  arma::cube  wi_store(nsave,NW,1);    wi_store.fill(0);
  arma::cube  vi_store(nsave,NV,1);    vi_store.fill(0);

  // STARTING VALUES
  arma::mat beta_draw   = beta_old;
  IntegerVector candind = seq_cpp(1,1000);


  // PROGRESS BAR
  Progress p(ntot, verbose==1);

  arma::mat y_w  = submat_arma(y_matrix, Wind);
  arma::mat y_v  = submat_arma(y_matrix, Vind);
  arma::mat Ni_w = submat_arma(Ni, Wind);
  arma::mat Ni_v = submat_arma(Ni, Vind);


  // GIBBS SAMPLER

for(int irep = 0; irep < ntot; irep++){

    // SAMPLE LATENT UTILITIES (LEMMA 2)

    arma::mat lambda   = exp(X * beta_draw);
    arma::mat lambda_w = submat_arma(lambda,Wind);
    arma::mat lambda_v = submat_arma(lambda,Vind);

    NumericVector Wi    = Rcpp::runif(NW, 0, 1);
    NumericVector Vi    = Rcpp::runif(NV, 0, 1);

    arma::mat exp_full   = (1/(Ni - y_matrix));
    arma::vec exponent2  = submat_arma(exp_full,Vind);
    arma::vec exponent1  = 1/submat_arma(y_matrix,Wind);

    NumericVector exp1 = as<NumericVector>(wrap(exponent1));
    NumericVector exp2 = as<NumericVector>(wrap(exponent2));

    arma::vec Wie = vecpow(Wi, exp1);
    arma::vec Vie = vecpow(Vi, exp2);

    arma::mat wi =   log((1 + lambda_w) % (1 / (Wie)) - lambda_w);
    arma::mat vi = - log(((1 + lambda_v) / lambda_v) % (1 / (Vie)) - (1/lambda_v));

    // SAMPLE SCALING FACTORS

      // SCALING FACTORS W

      arma::vec b_w = y_w + 1;
      arma::vec c_w = wi - log(lambda_w);

      NumericVector b_wn = as<NumericVector>(wrap(b_w));
      NumericVector c_wn = Rcpp::abs(as<NumericVector>(wrap(c_w)));

      arma::vec scaling_w = 1 / rcpp_pgdraw_b(b_wn, c_wn);


      // SCALING FACTORS V

      arma::vec b_v = Ni_v - y_v + 1;
      arma::vec c_v = vi - log(lambda_v);

      NumericVector b_vn  = as<NumericVector>(wrap(b_v));
      NumericVector c_vn  = Rcpp::abs(as<NumericVector>(wrap(c_v)));

      arma::vec scaling_v = 1 / rcpp_pgdraw_b(b_vn, c_vn);


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
    arma::vec wi_tilde = wi * sqrt(delta_star) + gamma_star;
    arma::vec vi_tilde = vi * sqrt(delta_star) + gamma_star;

    // DRAW GAMMA
    double gamma = 0;

    if(BOOST == 1){

      double u = R::runif(0,1);
      double L = 0;
      double O = 0;

      if(accu(y_matrix) == accu(Ni)){L = R_NegInf;}
      if(accu(y_matrix) != accu(Ni)){L = max(vi_tilde);}

      if(accu(y_matrix) == 0){O = R_PosInf;}
      if(accu(y_matrix) != 0){O = min(wi_tilde);}

      double O_ = O/sqrt(G0star);
      double L_ = L/sqrt(G0star);

      double prob  = u  * R::pt(O_, 2*d0, 1,0) + (1 - u) * R::pt(L_, 2*d0, 1, 0);
      gamma = sqrt(G0star) * R::qt(prob, 2*d0, 1, 0);

    }


    // DRAW COEFFICIENTS & DELTA

    arma::mat ma_w(N,1); ma_w.fill(0);
    arma::mat ma_v(N,1); ma_v.fill(0);
    arma::mat mb_w(N,1); mb_w.fill(0);
    arma::mat mb_v(N,1); mb_v.fill(0);
    arma::mat M_i(N,1);  M_i.fill(0);

    for(int ww = 0; ww < NW; ww++){
      ma_w.row(Wix(ww))  = wi_tilde(ww) / scaling_w(ww);
    }

    for(int vv = 0; vv < NV; vv++){
      ma_v.row(Vix(vv))  = vi_tilde(vv) / scaling_v(vv);
    }

    arma::mat ma_i     = ma_w + ma_v;
    arma::mat ma_mat   = X;

    for(int pp = 0; pp < P; pp++){
     ma_mat.col(pp) = ma_mat.col(pp) % ma_i;
    }

    arma::rowvec ma = arma::sum(ma_mat,0);


    for(int ww = 0; ww < NW; ww++){
        mb_w.row(Wix(ww))  = kappa_ww(ww);
    }

    for(int vv = 0; vv < NV; vv++){
        mb_v.row(Vix(vv))  = kappa_vv(vv);
    }

    arma::mat mb_i     = mb_w + mb_v;
    arma::mat mb_mat   = X;

    for(int pp = 0; pp < P; pp++){
      mb_mat.col(pp) = mb_mat.col(pp) % mb_i;
    }

    arma::rowvec mb = arma::sum(mb_mat,0);


    for(int ww = 0; ww < NW; ww++){
      M_i.row(Wix(ww))  = 1 / scaling_w(ww);
    }

    for(int vv = 0; vv < NV; vv++){
      M_i.row(Vix(vv))  = M_i.row(Vix(vv)) + 1/scaling_v(vv);
    }

    arma::mat X_star   = X;

    for(int pp = 0; pp < P; pp++){
      X_star.col(pp) = X_star.col(pp) % sqrt(M_i);
    }

    arma::mat Bn = inv(A0_inv + X_star.t() * X_star);


    // SAMPLE DELTA USING RESAMPLING TECHNIQUE

    double delta = 1;
    if(BOOST == 1){

      double dI     = d0 + 0.5 + 0.5 * (NV + NW);
      double SU     = sum(pow(wi_tilde,2) / scaling_w) + sum(pow(vi_tilde,2) / scaling_v);
      arma::vec aba = ma * Bn * ma.t();
      double maBnma = aba(0);
      arma::vec abb = ma * Bn * mb.t();
      double maBnmb = abb(0);

      double DI     = D0 + pow(gamma,2) / (2*G0) + 0.5 * SU - 0.5 * maBnma;
      double BI     = sum(wi_tilde % kappa_ww) + sum(vi_tilde % kappa_vv) - maBnmb;
      double dm     = d_m(DI,dI,BI);
      double Ip     = I_p(DI,dI,BI,dm);


      if(prior_ig == 1){

        double dI_star = -Ip * pow(dm,2) - 1;
        double DI_star = dm * (dI_star + 1);

        NumericVector delta_candidates = 1/Rcpp::rgamma(1000,dI_star, 1/DI_star);

        NumericVector diffpart  = (DI - DI_star) / delta_candidates;
        NumericVector diffpart2 = BI / sqrt(delta_candidates);
        NumericVector lweights = - (dI - dI_star) * log(delta_candidates) - diffpart + diffpart2;
        lweights = lweights - max(lweights);

        NumericVector weights = exp(lweights);
        NumericVector ws      = weights / sum(weights);

        int delta_index = sample(candind, 1, false, ws)[0];

        delta = delta_candidates(delta_index-1);

      }

      if(prior_ig == 0){

        double bI_star = -2 * Ip * pow(dm,2) - 1;
        double BI_star =  2 * (bI_star + 1) * sqrt(dm);

        NumericVector delta_candidates = 1/Rcpp::rgamma(1000,2*bI_star, 1/BI_star);
        delta_candidates = pow(delta_candidates,2);

        NumericVector dpart1 = DI / delta_candidates;
        NumericVector dpart2 = (BI+BI_star) / sqrt(delta_candidates);

        NumericVector lweights = - (dI - bI_star) * log(delta_candidates) - dpart1 + dpart2;

        lweights = lweights - max(lweights);
        NumericVector weights = exp(lweights);
        NumericVector ws      = weights / sum(weights);

        int delta_index = sample(candind, 1, false, ws)[0];

        delta = delta_candidates(delta_index-1);


      }

    }

     // SAMPLE COEFFICIENTS

     arma::mat m_i(N,1); m_i.fill(0);

     for(int ww = 0; ww < NW; ww++){
       m_i.row(Wix(ww))  = wi_tilde(ww) / scaling_w(ww) - sqrt(delta) * kappa_ww(ww);
     }

     for(int vv = 0; vv < NV; vv++){
       m_i.row(Vix(vv))  = m_i.row(Vix(vv)) + vi_tilde(vv) / scaling_v(vv) - sqrt(delta) * kappa_vv(vv);
     }


     arma::mat X_i = X;

     for(int pp = 0; pp < P; pp++){

       X_i.col(pp) = X_i.col(pp) % m_i;

     }

     arma::rowvec rX_i = arma::sum(X_i,0);
     arma::vec bn = Bn * rX_i.t();


       // DRAW FROM MVN

       arma::mat      sig =  delta * Bn;
       arma::vec      mu  =  bn;

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
         wi_store.row(irep - nburn)     = wi;
         vi_store.row(irep - nburn)     = vi;
         gamma_store.row(irep - nburn)  = gamma;
         delta_store.row(irep - nburn)  = delta;

        }

    // PROGRESS BAR

    if(verbose == 1) p.increment(); // update progress

    // CHECK USER INTERRUPTION

    if(irep % 20 == 0) Rcpp::checkUserInterrupt();

 }

return List::create(Named("beta")     = beta_store,
                    Named("v.i")      = vi_store,
                    Named("w.i")      = wi_store,
                    Named("gamma")    = gamma_store,
                    Named("delta")    = delta_store);


}











