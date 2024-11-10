upg.mnl   =   function(y.matrix,
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

  # PRELIMINARIES

  y.categ  = max.col(y.matrix)

  # DIMENSIONS

  N        = nrow(X)
  P        = ncol(X)
  K        = ncol(y.matrix)

  # GIBBS SAMPLER PARAMETERS

  nsave    = nsave
  nburn    = nburn
  ntot     = nburn + nsave

  # PRIOR ON COEFFICIENTS AND INTERCEPT

  A0.inv   = diag(1 / c(A0, rep(B0, P-1)), P)

  # STORAGE

  beta.store     = array(NA, c(nsave, P, K - 1))

  # STARTING VALUES

  if(is.null(beta.start)) {beta.draw = array(0, c(P, K-1))} else {beta.draw = beta.start[, -K]}
  z = matrix(0, N, K)

  # VERBOSE

  if(verbose) pb = txtProgressBar(min = 0, max = ntot, style = 3)

  # LOG INVERSE GAMMA DENSITY
  dinvgamma.log = function(x, a, b){

    a * log(b) - lgamma(a) + (a + 1) * log(1 / x) - b/x

  }

  # - - - - START MCMC ALGORITHM

  for(irep in 1:ntot){

    for(kk in 1:(K-1)){ # - START LOOP OVER ALL CATEGORIES BUT BASELINE


      # - - - - STEP 1: SAMPLE LATENT UTILITIES

      lambda      = cbind(exp(X %*% beta.draw), 1)
      lambda.star = rowSums(lambda)
      Ui          = runif(N)
      Vi          = matrix(runif(N * K), N, K)
      ut.matrix   = - log(Ui) / lambda.star - (log(Vi) / lambda) * (1 - y.matrix)
      y.d         = - log(ut.matrix)
      z           = y.d[,kk] - matrixStats::rowMaxs(y.d[,-kk,drop=F])
      xi          = log(rowSums(lambda[,-kk,drop=F]))

      # - - - - STEP 2: SAMPLE LATENT SCALING FACTORS

      omega       = pgdraw::pgdraw(2, abs(z + xi - log(lambda[,kk])))


      # - - - - STEP 3: SHIFT MOVE

      # DRAW FROM WORKING PRIOR

      gamma.star = rnorm(1, 0, sqrt(G0))

      if(!gamma.boost){gamma.star   = 0}

      # SHIFT UTILITIES

      ztilde     = z + gamma.star

      # POSTERIOR QUANTITIES FOR GAMMA CONDITIONAL
      XW         = X * omega
      tXW        = t(XW)
      Bn         = chol2inv(chol(A0.inv + tXW %*% X))

      if(gamma.boost){

        mi         = (ztilde + xi) * omega# - 0.5
        mb         = colSums(XW)
        mg         = sum(mi)
        mn         = t(X) %*% mi
        tmbBn      = t(mb) %*% Bn
        Gn         = 1 / ( (1 / G0) + sum(omega) - tmbBn %*% mb)
        gn         = Gn * (mg - tmbBn %*% mn)

        # TRUNCATION POINTS FOR GAMMA CONDITIONAL

        U          = ifelse(sum(y.matrix[,kk]) == 0,  Inf,  min(ztilde[y.matrix[, kk] == 1]))
        L          = ifelse(sum(y.matrix[,kk]) == N, -Inf,  max(ztilde[y.matrix[, kk] == 0]))

        # SIMULATE FROM GAMMA POSTERIOR

        gamma      = truncnorm::rtruncnorm(1, a = L, b = U, mean = gn, sd = sqrt(Gn))

      }

      if(!gamma.boost){gamma=0}

      # REVERSE SHIFT
      z.shift    = ztilde - gamma



      # - - - - STEP 4: SCALE MOVE

      # DRAW FROM WORKING PRIOR
      if(!delta.boost){delta.star  = 1}

      # DRAW FROM DELTA CONDITIONAL USING RESAMPLING METHOD BASED ON IG AUXILIARY PRIOR
      if(delta.boost){

        delta.star = 1/rgamma(1, d0, D0)


        # POSTERIOR QUANTITIES COEFFICIENTS
        ma          = t(X) %*% ((z.shift) * omega)
        mb          = t(X) %*% (omega * (xi))

        dI          = c(d0 + 0.5 * N)
        DI          = c(D0 + 0.5 * delta.star * (sum((z.shift)^2 * omega) - t(ma) %*% Bn %*% ma))
        BI          = - c(sqrt(delta.star) * (sum(z.shift * omega * (xi)) - t(ma) %*% Bn %*% mb))

        # MODE AND CURVATURE OF POSTERIOR
        dm = (16 * (DI^2)) / (BI + sqrt(BI^2 + 16 * DI * (dI +1)))^2
        IP = - sqrt(BI^2 + 16 * DI * (dI + 1)) / (4 * (dm)^(5/2))

        # MOMENTS OF AUXILIARY PRIOR
        dI.star = - IP * (dm^2) - 1
        DI.star = dm * (dI.star + 1)

        # DRAW CANDIDATES FROM MATCHED AUXILIARY PRIOR (INVERSE GAMMA)

        delta.candidates = 1 / rgamma(30, dI.star, DI.star)
        lprior           = dinvgamma.log(delta.candidates, dI.star, DI.star)

        # RESAMPLING

        lweights         = - (dI + 1) * log(delta.candidates) - DI / delta.candidates + BI / sqrt(delta.candidates) - lprior
        lweights         = lweights - max(lweights)
        weights          = exp(lweights)
        weights          = weights / sum(weights)

        # DRAW DELTA

        delta.index = sample(1:length(delta.candidates), 1, prob = weights)
        delta       = delta.candidates[delta.index]


      }

      if(!delta.boost){delta      = 1}

      # - - - - STEP 5: SAMPLE COEFFICIENTS

      mb               = t(X) %*% (omega * (xi))
      ma               = t(X) %*% (omega * (sqrt(delta.star / delta) * z.shift))
      mk               = mb + ma

      bn               = (Bn %*% mk)
      sqrtBn           = t(chol(Bn))
      beta.draw[,kk]   = bn + sqrtBn %*% rnorm(P)

    } # - END LOOP OVER ALL CATEGORIES BUT BASELINE


    # - STORE POSTERIOR SAMPLES AFTER BURN-IN PERIOD

    if(irep>nburn){

      beta.store[irep-nburn,,]            = beta.draw

    }

    if(verbose==T) {setTxtProgressBar(pb, irep)}

  }

  # - - - - END MCMC ALGORITHM

  # ADD BASELINE CATEGORY (RESTRICTED TO ZERO)

  beta.store.bl       = array(0, c(nsave, P, K))
  beta.store.bl[,,-K] = beta.store
  beta.store          = beta.store.bl

  return(list(y = y.matrix, X = X, beta = beta.store))

}
