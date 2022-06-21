upg.probit = function(y,
                      X,
                      nsave       = 10000,
                      nburn       = 2000,
                      d0          = 2.5,
                      D0          = 1.5,
                      G0          = 100,
                      B0          = 4,
                      A0          = 4,
                      gamma.boost = T,
                      delta.boost = T,
                      beta.start  = NULL,
                      verbose     = T)
{

  # - - - - MCMC SETUP

  # DIMENSIONS

  N        = nrow(X)
  P        = ncol(X)

  # GIBBS SAMPLER PARAMETERS

  nsave    = nsave
  nburn    = nburn
  ntot     = nburn + nsave

  # PRIOR ON COEFFICIENTS AND INTERCEPT

  A0.inv   = diag(1 / c(A0, rep(B0, P-1)), P)

  # STORAGE

  beta.store     = array(NA, c(nsave, P))
  gamma.store    = array(NA, c(nsave, 1))
  delta.store    = array(NA, c(nsave, 1))

  # STARTING VALUES

  if(is.null(beta.start)) {beta.draw = array(0, c(P, 1))} else {beta.draw = beta.start}

  # FUNCTIONS FOR SAMPLING LATENT UTILITIES

  Fe_ = function(p){qnorm(p)}
  Fe  = function(x){pnorm(x)}

  # VERBOSE

  if(verbose) pb = txtProgressBar(min = 0, max = ntot, style = 3)

  # SOME PRECOMPUTATIONS

  Sy0   = sum(y) == 0
  SyN   = sum(y) == N
  y0    = y == 0
  y1    = y == 1
  omega = rep(1, N)

  # - - - - START MCMC ALGORITHM

  for(irep in 1:ntot){

    # - STEP 1: SAMPLE LATENT UTILITIES

    U.i    = runif(N)
    lambda = X %*% beta.draw
    z      = lambda + Fe_(y + U.i * (1 - y - (1 - Fe( -lambda))))


    # - STEP 2: SAMPLE LATENT SCALING FACTORS (ALL EQUAL 1 IN PROBIT)


    # - STEP 3: SHIFT MOVE

    # DRAW FROM WORKING PRIOR

    gamma.star = rnorm(1, 0, sqrt(G0))

    if(!gamma.boost){gamma.star   = 0}

    # SHIFT UTILITIES

    ztilde     = z + gamma.star

    # POSTERIOR QUANTITIES FOR GAMMA CONDITIONAL

    XW         = X * omega
    tXW        = t(XW)
    Bn         = chol2inv(chol(A0.inv + tXW %*% X)) # could be precomputed for probit
    mb         = colSums(XW)
    mg         = sum(ztilde * omega)
    mn         = tXW %*% ztilde
    tmbBn      = t(mb) %*% Bn
    Gn         = 1 / ( (1 / G0) + N - tmbBn %*% mb)
    gn         = Gn * (mg - tmbBn %*% mn)

    # TRUNCATION POINTS FOR GAMMA CONDITIONAL

    U          = ifelse(Sy0,  Inf,  min(ztilde[y1]))
    L          = ifelse(SyN, -Inf,  max(ztilde[y0]))

    # SIMULATE FROM GAMMA POSTERIOR

    gamma      = truncnorm::rtruncnorm(1, a = L, b = U, mean = gn, sd = sqrt(Gn))

    if(!gamma.boost){gamma=0}

    # REVERSE SHIFT

    z.shift    = ztilde - gamma


    # - STEP 4: SCALE MOVE

    # DRAW FROM WORKING PRIOR

    delta.star = 1/rgamma(1, d0, D0)

    if(!delta.boost){delta.star  = 1}

    # POSTERIOR QUANTITIES FOR DELTA CONDITIONAL

    mn = tXW %*% z.shift
    bn = Bn %*% mn

    # SIMULATE FROM DELTA POSTERIOR

    delta      = 1 / rgamma(1, d0 + 0.5 * N, D0 + 0.5 * delta.star * (t(bn) %*% A0.inv %*% bn + sum(omega * (z.shift - X %*% bn)^2)))

    if(!delta.boost){delta      = 1}


    # - STEP 5: SAMPLE COEFFICIENTS

    sqrtBn    = t(chol(Bn))
    beta.draw = sqrt(delta.star / delta) * bn + sqrtBn %*% rnorm(P)


    # - STORE POSTERIOR SAMPLES AFTER BURN-IN PERIOD

    if(irep>nburn){

      beta.store[irep-nburn,]            = beta.draw
      gamma.store[irep-nburn,]           = gamma.star - gamma
      delta.store[irep-nburn,]           = delta.star / delta

    }

    if(verbose==T) {setTxtProgressBar(pb, irep)}

  }

  # - - - - END MCMC ALGORITHM


  return(list(y = y, X = X, beta = beta.store, gamma = gamma.store, delta = delta.store))

}
