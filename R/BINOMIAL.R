upg.binomial = function(y        = NULL,
                        X        = NULL,
                        Ni       = NULL,
                        BOOST    = T,
                        verbose  = F,
                        nburn    = 1000,
                        nsave    = 1000,
                        A0       = 1,
                        d0       = 0.5,
                        D0       = 0.5,
                        G0       = 99,
                        beta.old = NULL,
                        prior.ig = T){



  ############################
  # - GIBBS SAMPLER SETUP -  #
  ############################

  if(is.null(beta.old)) {beta.draw = array(0, c(ncol(X),1))} else {beta.draw = beta.old}
  BOOST    = ifelse(BOOST, 1,0)
  verbose  = ifelse(verbose, 1,0)
  prior_ig = ifelse(prior.ig, 1, 0)

  ############################
  # - GIBBS SAMPLING      -  #
  ############################

  res =  sample_binom(y_matrix = y,
                      Ni       = Ni,
                      X        = X,
                      beta_old = beta.draw,
                      nsave    = nsave,
                      nburn    = nburn,
                      A0       = A0,
                      d0       = d0,
                      D0       = D0,
                      G0       = G0,
                      BOOST    = BOOST,
                      verbose  = verbose,
                      prior_ig = prior_ig)

  #STORAGE
  beta.store     = res$beta
  gamma.store    = res$gamma
  delta.store    = res$delta
  wi.store       = res$w.i[,,1]
  vi.store       = res$v.i[,,1]

  return(list(beta        = beta.store,
              w.i         = wi.store,
              v.i         = vi.store,
              gamma       = gamma.store,
              delta       = delta.store))

}
