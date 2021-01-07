upg.mnl    = function(y        = NULL,
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

  if(0 %in% y){y = y+1}
  if(length(unique(y)==1)) {y.matrix = matrix(0, nrow=nrow(y), ncol = 2)}
  if(length(unique(y))>1)  {y.matrix = matrix(0, nrow=nrow(y), ncol = length(unique(y)))}
  for(ii in 1:nrow(y)) y.matrix[ii,y[ii,]] = 1
  y.categ  = as.matrix(max.col(y.matrix))-1
  if(is.null(beta.old)) {beta.draw = array(0, c(ncol(X),ncol(y.matrix)))} else {beta.draw = beta.old}
  BOOST   = ifelse(BOOST, 1,0)
  verbose = ifelse(verbose, 1,0)

  ############################
  # - GIBBS SAMPLING      -  #
  ############################

res =  sample_mnl(y_matrix = y.matrix,
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
  gamma.store    = array(res$gamma,c(nrow(res$gamma),1,ncol(res$gamma)))
  delta.store    = array(res$delta,c(nrow(res$delta),1,ncol(res$delta)))
  y.latent.store = res$y.latent


  return(list(beta        = beta.store,
              y.latent    = y.latent.store,
              gamma       = gamma.store,
              delta       = delta.store))

}
