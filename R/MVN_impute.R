MVN_impute  <- function( Y, reltol=1e-4, intercept=TRUE, maxit=1e2, trace=FALSE ){
  
  if( any( colMeans( is.na(Y) ) == 1 ) )
    stop( 'Completely blank columns of Y are not allowed' )

  N <- nrow(Y)
  P <- ncol(Y)    
  miss.indices	<- sapply( 1:N, function(n) which( is.na(Y[n,]) ) )
  nmiss		<- sapply( miss.indices, length )
  
  if( sum(nmiss) == 0 ) 
    return(list( Y=Y, Sigma=var(Y) ))
  
  loopers   <- which( nmiss > 0 & nmiss < P )
  emptyInds <- which( nmiss==P )
  fullInds  <- which( nmiss==0 )
  nempty    <- length(emptyInds)
  
  ######## initialization
  if( intercept ){
    mu		= colMeans( Y, na.rm=TRUE )
  } else {
    mu  <- rep( 0, P )
  }
  
  Yhat	= Y
  for(p in 1:P)	##### initialize missing values to column means
    Yhat[which( is.na(Yhat[,p])),p] = mu[p]
  
  Sigma   = .9*var(Yhat)+.1*diag(P) # identity causes weird singularity and var(Yhat) needn't be pd
  
  if( trace ){
    cat( sprintf("%10s, %12s\n",  "counter", '% imp change') )
    cat( sprintf("%10d, %12.2e\n", 0,         NA) )
  }

  for( it in 1:maxit ){
    Y.old   <- Yhat
    
    ####### E step
    Yhat		<- Y
    Sig.n	<- array( 0, dim=c(N,P,P) ) ## for M step
    for( n in loopers ){
      m   	    <- miss.indices[[n]]
      o   	    <- (1:P)[-m]	    
      Yhat[n,m]	    <- mu[m] + Sigma[m,o,drop=FALSE]%*%solve(Sigma[o,o])%*%t(Y[n,o,drop=FALSE]-mu[o])
      Sig.n[n,m,m]  <- schur( Sigma, m ) ### used in M
    }
    
    ####### E step cont'd
    if( nempty > 0 ){
      Yhat[emptyInds,]	<- matrix( mu, nempty, P, byrow=TRUE )
      for( n in emptyInds )
        Sig.n[n,,]	<- Sigma ### used in M
    }

    ##### determine convergence
    delta <- 100 * mean(( (Yhat-Y.old)[is.na(Y)] )^2 ) / mean( ( Yhat[is.na(Y)] )^2 )
    if( trace )
      cat( sprintf("%10d, %12.2e\n", it, delta ) )
    if ( delta/100 < reltol )
      break
    
    ####### M step
    if( intercept )
      mu	<- colMeans( Yhat )
    
    mu.mat    <- matrix( mu, N, P, byrow=TRUE )
    S.part.1  <- 1/N * t(Yhat - mu.mat) %*% (Yhat - mu.mat)
    Sigma     <- S.part.1 + 1/N * apply( Sig.n, 2:3, sum )
  }
  return(list(Y=Yhat,Sigma=Sigma))
}
