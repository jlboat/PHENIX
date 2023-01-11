lmm_obj <- function(delta,yprime,Lam.K){
  denominator <- Lam.K+delta
  sig_g       <- mean( yprime^2 / denominator )
  sum(dnorm(yprime,mean=0,sd=sqrt(sig_g*(denominator)),log=TRUE))
}

fit_lmm <- function(y,yprime,K,eig.K,Lam.K,constrained=TRUE){

  if( missing( yprime ) | missing( Lam.K ) ){
    miss  <- which(is.na(y))
    if( length( miss ) > 0 ){
      y <- y[ -miss ]
      K <- K[ -miss, -miss ]
    }
    if( missing( eig.K ) ){
      eig.K <- eigen( K, symmetric=TRUE )
    } else {
      if( nrow(eig.K$vec) != nrow(K) )
        stop(paste0( 'eig.K has ', nrow(eig.K$vec), ' samples but # non-missing samples=', length( y ) ))
    }
    yprime  <- t( eig.K$vec ) %*% y
    Lam.K   <- eig.K$val
    rm( eig.K, y )
  }

  if( constrained ){
    opt_lims  <- c(.001,1000)
    delta     <- optimise(lmm_obj,opt_lims,maximum=TRUE,Lam.K=Lam.K,yprime=yprime)$maximum
  } else {

    opt_lims_plus <- c( -min(Lam.K), 1000 )
    opt_lims_minus<- c( -1000, -max(Lam.K) )

    opt_plus  <- optimise(lmm_obj,opt_lims_plus,  maximum=TRUE,Lam.K=Lam.K,yprime=yprime)
    opt_minus <- optimise(lmm_obj,opt_lims_minus, maximum=TRUE,Lam.K=Lam.K,yprime=yprime)

    delta     <- ifelse(
      opt_plus$objective > opt_minus$objective,
      opt_plus$maximum,
      opt_minus$maximum
    )

  }

  sig_g       <- mean( yprime^2/(Lam.K+delta) )
  sig_e       <- delta*sig_g
  return( list(sig_g=sig_g,sig_e=sig_e,h2=sig_g/(sig_g+sig_e)) )
}

lmm.impute <- function( Y, K, eig.marg, trace=FALSE ) {
  Yhat    <- Y
  U       <- Y + NA
  h2      <- rep( NA, ncol(Y) )
  for( p in 1:ncol(Y) ){
  
    if( trace )
      print( paste0( 'Running LMM for phenotype #', p ) )
  
    m <- which( is.na(Y[,p]) )
    o <- which( ! is.na(Y[,p]) )

    if( length( o ) == 0 ){
      warning(paste0( 'Phenotype ', p, ' is completely blank; imputing to 0' ))
      h2  [p]   <- NA
      Yhat[,p]  <- 0
      U   [,p]  <- 0
      next 
    }

    if( missing(eig.marg) ){
      out <- fit_lmm(y=Y[,p],K=K)
    } else {
      out <- fit_lmm( yprime=t( eig.marg[[p]]$vec ) %*% Y[o,p], Lam.K=eig.marg[[p]]$val )
    }
    sig_g <- out$sig_g
    sig_e <- out$sig_e
    h2[p] <- sig_g / ( sig_g+sig_e )
    
    ### Compute breeding values/predicted phenotypes
    Sigma_mo  <- sig_g*K[,o,drop=FALSE]
    Sigma_oo  <- solve( sig_g*K[o,o,drop=FALSE] + sig_e*diag(length(o)) )
    U[,p]     <- Sigma_mo %*% Sigma_oo %*% Y[o,p,drop=FALSE]
    if( length(o) < nrow(Y) )
      Yhat[-o,p] <-U[-o,p]

  }
  return( list( Yhat=Yhat, h2=h2, U=U ) )
}
