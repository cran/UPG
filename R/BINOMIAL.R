upg.binomial = function(y,
                        X,
                        Ni,
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

  # THINGS FOR SAMPLING LATENT UTILITIES

  # SOME INDICATORS
  W.ind = y > 0
  V.ind = y < Ni
  NW    = sum(W.ind)
  NV    = sum(V.ind)

  # SOME PRECOMPUTATIONS
  kappa.w = (1 - y) / 2
  kappa.v = (Ni - y - 1) / 2


  # INVERSE GAMMA LOG DENSITY
  log.IG = function(x, a, b){a * log(b) - lgamma(a) + (a + 1) * log(1 / x) - b/x}

  # VERBOSE

  if(verbose) pb = txtProgressBar(min = 0, max = ntot, style = 3)

  post.delta = function(delta, dI, DI, BI){

    (dI+1) * log(1/delta) - DI/delta + (BI / sqrt(delta))

  }

  # - - - - START MCMC ALGORITHM

  for(irep in 1:ntot){

    # - STEP 1: SAMPLE LATENT UTILITIES USING THEOREM 3

    lambda = exp(X %*% beta.draw)
    Wi     = runif(NW)
    Vi     = runif(NV)
    wi     =   log((1 + lambda[W.ind]) * (1 / (Wi^(1/y[W.ind]))) - lambda[W.ind])
    vi     = - log(((1 + lambda[V.ind]) / lambda[V.ind]) * (1 / (Vi^(1/(Ni - y)[V.ind]))) - (1/lambda[V.ind]))


    # - STEP 2: SAMPLE LATENT SCALING FACTORS

    omega.w = pgdraw::pgdraw(y[W.ind] + 1, abs(wi - log(lambda[W.ind])))
    omega.v = pgdraw::pgdraw(Ni[V.ind] - y[V.ind] + 1, abs(vi - log(lambda[V.ind])))


    # - STEP 3: SHIFT MOVE

    # DRAW FROM WORKING PRIOR

    gamma.star = rnorm(1, 0, sqrt(G0))

    if(!gamma.boost){gamma.star   = 0}

    # SHIFT UTILITIES

    wtilde     = wi + gamma.star
    vtilde     = vi + gamma.star

    # POSTERIOR QUANTITIES FOR GAMMA CONDITIONAL

    Mi         = numeric(length = N)
    Mi[W.ind]  = omega.w
    Mi[V.ind]  = omega.v + Mi[V.ind]

    XW         = X * Mi
    tXW        = t(XW)
    Bn         = chol2inv(chol(A0.inv + tXW %*% X))

    mi         = numeric(length = N)
    mi[W.ind]  = wtilde * omega.w - kappa.w[W.ind]
    mi[V.ind]  = vtilde * omega.v - kappa.v[V.ind] + mi[V.ind]

    mb         = colSums(XW)
    mg         = sum(mi)
    mn         = t(X) %*% mi
    tmbBn      = t(mb) %*% Bn
    Gn         = 1 / ( (1 / G0) + sum(Mi) - tmbBn %*% mb)
    gn         = Gn * (mg - tmbBn %*% mn)

    # TRUNCATION POINTS FOR GAMMA CONDITIONAL

    U          = min(wtilde)
    L          = max(vtilde)

    # SIMULATE FROM GAMMA POSTERIOR

    gamma      = truncnorm::rtruncnorm(1, a = L, b = U, mean = gn, sd = sqrt(Gn))

    if(!gamma.boost){gamma=0}

    # REVERSE SHIFT

    w.shift    = wtilde - gamma
    v.shift    = vtilde - gamma


    # - STEP 4: SCALE MOVE

    # DRAW FROM WORKING PRIOR

    delta.star = 1/rgamma(1, d0, D0)

    if(!delta.boost){delta.star  = 1}

    # DRAW FROM DELTA CONDITIONAL USING RESAMPLING METHOD BASED ON IG AUXILIARY PRIOR

    # POSTERIOR QUANTITIES COEFFICIENTS

    ma.y        = numeric(length = N)
    ma.y[W.ind] = w.shift * omega.w
    ma.y[V.ind] = v.shift * omega.v + ma.y[V.ind]
    ma          = t(X) %*% ma.y

    mb.y        = numeric(length = N)
    mb.y[W.ind] = kappa.w[W.ind]
    mb.y[V.ind] = kappa.v[V.ind] + mb.y[V.ind]
    mb          = t(X) %*% mb.y

    # COMPUTE RELEVANT LIKELIHOOD QUANTITIES

    dI = c(d0 + 0.5 * (NV + NW))
    DI = c(D0 + 0.5 * delta.star * ((sum((w.shift)^2 * omega.w) + sum((v.shift)^2 * omega.v)) - t(ma) %*% Bn %*% ma))
    BI = c(sqrt(delta.star) * (sum(w.shift * kappa.w[W.ind]) + sum(v.shift * kappa.v[V.ind]) - t(ma) %*% Bn %*% mb))


    # MODE AND CURVATURE OF POSTERIOR

    dm = (16 * (DI^2)) / (BI + sqrt(BI^2 + 16 * DI * (dI +1)))^2
    IP = - sqrt(BI^2 + 16 * DI * (dI + 1)) / (4 * (dm)^(5/2))

    # MOMENTS OF AUXILIARY PRIOR

    dI.star = - IP * (dm^2) - 1
    DI.star = dm * (dI.star + 1)

    # x1 = seq(0,10,0.01)
    # y1 = post.delta(x1, dI, DI, BI)
    # y2 = dgamma(x1, dI.star, DI.star,log=T)
    #
    # plot(x1, (y1), type = "l")
    # lines(x1, y2, type="l", col="blue")
    #   abline(v=dm)
    # abline(v = (DI.star / (dI.star + 1)), col="red")
    # abline(v = (DI / (dI.star + 1)), col="blue")
    # abline(v = (DI.star / (dI + 1)), col="green")

    # RESAMPLING

    delta.candidates = 1 / rgamma(100, dI.star, DI.star)
    #delta.candidates = runif(1000, 0, 10)
    lprior   = log.IG(delta.candidates, dI.star, DI.star)
    #lprior    = dunif(delta.candidates,0,10,log=T)
    lweights = - (dI + 1) * log(delta.candidates) - DI / delta.candidates + BI / sqrt(delta.candidates) - lprior
    #  lweights = - (dI - dI.star) * log(delta.candidates) - (DI - DI.star) / delta.candidates + BI / sqrt(delta.candidates)
    lweights = lweights - max(lweights)
    weights  = exp(lweights)
    weights  = weights / sum(weights)

    # DRAW DELTA

    delta.index = sample(1:100, 1, prob = weights)
    delta       = delta.candidates[delta.index]
    if(!delta.boost){delta      = 1}

    # - STEP 5: SAMPLE COEFFICIENTS
    bn          = Bn %*% ((sqrt(delta.star / delta) * ma - mb))
    sqrtBn      = t(chol(Bn))
    beta.draw   = bn + sqrtBn %*% rnorm(P)


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
