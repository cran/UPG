upg.logit  = function(y        = NULL,
                      X        = NULL,
                      nsave    = 500,
                      nburn    = 500,
                      A0       = 1,
                      d0       = 0.5,
                      D0       = 0.5,
                      G0       = 99,
                      BOOST    = T,
                      verbose  = F,
                      beta.old = NULL)
{




  ############################
  # - GIBBS SAMPLER SETUP -  #
  ############################

  y.matrix = matrix(0, nrow=nrow(y), ncol=2)
  for(ii in 1:nrow(y)) y.matrix[ii,ifelse(y[ii,]==1,1,2)] = 1
  y.categ  = as.matrix(max.col(y.matrix))
  if(is.null(beta.old)) {beta.draw = array(0, c(ncol(X),1))} else {beta.draw = beta.old}
  BOOST   = ifelse(BOOST, 1,0)
  verbose = ifelse(verbose, 1,0)

  ############################
  # - GIBBS SAMPLING      -  #
  ############################

res =  sample_logit(y_matrix = y.matrix,
                y_categ  = y.categ,
                X        = X,
                beta_old = beta.draw,
                nsave    = nsave,
                nburn    = nburn,
                A0       = A0,
                d0       = d0,
                D0       = D0,
                G0       = G0,
                BOOST    = BOOST,
                verbose  = verbose)

  #STORAGE
  beta.store     = res$beta
  gamma.store    = array(res$gamma,c(nrow(res$gamma), ncol(res$gamma),1))
  delta.store    = array(res$delta,c(nrow(res$delta), ncol(res$delta),1))
  y.latent.store = res$y.latent


  return(list(beta        = beta.store,
              y.latent    = y.latent.store,
              gamma       = gamma.store,
              delta       = delta.store))

}
