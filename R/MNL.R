upg.mnl   = function(y.matrix,
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

  # FUNCTION FOR SAMPLING LATENT UTILITIES

  sample.util.new = function(lambda.star,
                             lambda.k,
                             lambda.a,
                             y.matrix,
                             categ,
                             y.categ,
                             K){


    N   = length(lambda.star)
    Ui  = runif(N)
    Vki = runif(N)
    V0i = runif(N)
    Vai = runif(N)

    ut.matrix   = matrix(nrow = N, ncol = 3)

    #first part of difference
    utilities = - log(Ui) / lambda.star

    #individual is categ
    ut.matrix[,1] = utilities - (log(Vki) / lambda.k) * as.numeric(y.categ != categ)

    #individual is other
    ut.matrix[,2] = utilities - (log(Vai) / lambda.a) * as.numeric((y.categ %in% c(categ,K)))

    #individual is baseline
    ut.matrix[,3] = utilities - log(V0i) * as.numeric(y.categ != K)

    return(-log(ut.matrix))

  }

  # VERBOSE

  if(verbose) pb = txtProgressBar(min = 0, max = ntot, style = 3)

  # - - - - START MCMC ALGORITHM

  for(irep in 1:ntot){

    for(kk in 1:(K-1)){ # - START LOOP OVER ALL CATEGORIES BUT BASELINE


      # - STEP 1: SAMPLE LATENT UTILITIES

      lambda      = cbind(exp(X %*% beta.draw), 1)
      lambda.star = rowSums(lambda)
      lambda.k    = lambda[,kk]
      lambda.a    = rowSums(lambda[,-c(kk, ncol(lambda)), drop=F])

      y.d         = sample.util.new(lambda.star,
                                    lambda.k,
                                    lambda.a,
                                    y.matrix = y.matrix,
                                    categ    = kk,
                                    K        = K,
                                    y.categ  = y.categ)

      # DIFFERENCE OF ACTIVE CATEGORY AND BASELINE

      z         = y.d[, 1] - y.d[, 3]

      # DIFFERENCE OF ACTIVE CATEGORY TO MAXIMUM

      allothers = y.d[, -1]
      diffUtil  = y.d[,  1] - allothers[cbind(seq_len(nrow(allothers)), max.col(allothers))]


      # - STEP 2: SAMPLE LATENT SCALING FACTORS

      omega  = pgdraw::pgdraw(2, abs(z - log(lambda.k)))


      # - STEP 3: SHIFT MOVE

      # DRAW FROM WORKING PRIOR

      gamma.star = rnorm(1, 0, sqrt(G0))

      if(!gamma.boost){gamma.star   = 0}

      # SHIFT UTILITIES

      ztilde     = z        + gamma.star
      diffUtil   = diffUtil + gamma.star

      # POSTERIOR QUANTITIES FOR GAMMA CONDITIONAL

      XW         = X * omega
      tXW        = t(XW)
      Bn         = chol2inv(chol(A0.inv + tXW %*% X))
      mb         = colSums(XW)
      mg         = sum(ztilde * omega)
      mn         = tXW %*% ztilde
      tmbBn      = t(mb) %*% Bn
      Gn         = 1 / ( (1 / G0) + sum(omega) - tmbBn %*% mb)
      gn         = Gn * (mg - tmbBn %*% mn)

      # TRUNCATION POINTS FOR GAMMA CONDITIONAL

      U          = ifelse(sum(y.matrix[,kk]) == 0,  Inf,  min(diffUtil[y.matrix[, kk] == 1]))
      L          = ifelse(sum(y.matrix[,kk]) == N, -Inf,  max(diffUtil[y.matrix[, kk] == 0]))

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

      sqrtBn          = t(chol(Bn))
      beta.draw[, kk] = sqrt(delta.star / delta) * bn + sqrtBn %*% rnorm(P)


    } # - END LOOP OVER ALL CATEGORIES BUT BASELINE


    # - STORE POSTERIOR SAMPLES AFTER BURN-IN PERIOD

    if(irep>nburn){

      beta.store[irep-nburn,,]            = beta.draw

    }

    if(verbose==T) {setTxtProgressBar(pb, irep)}

  }

  # - - - - END MCMC ALGORITHM

  # add zero dimension for baseline coefficients

  beta.store.bl       = array(0, c(nsave, P, K))
  beta.store.bl[,,-K] = beta.store
  beta.store          = beta.store.bl



  return(list(y = y.matrix, X = X, beta = beta.store))

}
