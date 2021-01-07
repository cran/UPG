#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace Rcpp;



arma::vec Fe_(arma::vec x){

  NumericVector x1 = as<NumericVector>(wrap(x));
  arma::vec qnormed = Rcpp::qnorm(x1, 0, 1, 1, 0);
  return qnormed;

  }

arma::vec Fe(arma::vec x){

  NumericVector x1 = as<NumericVector>(wrap(x));
  arma::vec pnormed = Rcpp::pnorm(x1, 0, 1, 1, 0);
  return pnormed;

  }


//' @name sample_probit
//' @noRd
// [[Rcpp::export]]
List sample_probit(arma::mat y_matrix,
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


  // CREATE ALL NECESSARY THINGS
  arma::mat U_i(N,1);
  arma::mat lambda(N,1);
  arma::mat y_d(N,1);

  // GIBBS SAMPLER

  for(int irep = 0; irep < ntot; irep++){

    // SAMPLE LATENT UTILITIES

    arma::mat U_i    = Rcpp::runif(N, 0, 1);
    arma::mat lambda = X * beta_draw;
    arma::mat y_d    = lambda + Fe_(y_matrix.col(0) + U_i % (1 - y_matrix.col(0) - (1 - Fe(-lambda))));


    // SAMPLE SCALING FACTORS

    arma::vec omega(N); omega.fill(1);

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



