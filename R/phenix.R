phenix <- function ( Y, K, Q, lam_K, test=FALSE, test.frac=.05, test.cols=1:ncol(Y), seed=8473, quantnorm=FALSE, scale=TRUE, trim=FALSE, trim.sds=4, ... ) 
{

  set.seed( seed )

  intercept <- all( Y[,1] == 1 )
  if( test ){
    t0      <- proc.time()[3]
    obs     <- which( !is.na(Y) )
    if( any( test.cols < 1 ) | any( test.cols > ncol(Y) ) )
      stop( 'test.cols must be an integer vector with entries in 1, 2, ..., ncol(Y)' )
    if( test.frac <= 0 | test.frac >= 1 )
      stop( 'test.frac must be between 0 and 1' )
    if( intercept & ( 1 %in% test.cols ) ){
      if(!all( test.cols == 1:ncol(Y) ))
        warning( 'Testing intercept not allowed' )
      test.cols <- test.cols[ - which( test.cols == 1 ) ]
    }
    obs     <- intersect( obs, c( matrix( 1:(nrow(Y)*ncol(Y)), nrow(Y), ncol(Y) )[ ,test.cols ] ) )
    mask    <- sample( obs, floor( test.frac*length( obs ) ) )
    Y.input <- Y
    Y[mask] <- NA
  }

  if( trim ){
    outlier.mat <- matrix( 0, nrow(Y), ncol(Y) )
    for( p in 1:ncol(Y) ){
      y         <- Y[,p]
      mu        <- mean(y,na.rm=TRUE)
      sig       <- sd(y,na.rm=TRUE)
      if( sig == 0 & p == 1 )
        next
      outliers_p  <- which( y < mu - trim.sds*sig | y > mu + trim.sds*sig )
      if( length( outliers_p ) > 0 )
        outlier.mat[ outliers_p, p ]  <- 1 
      rm( outliers_p )
    }
    Y_outliers  <- Y
    outliers    <- which( outlier.mat == 1 )
    if( length( outliers ) > 0 )
      Y[ outliers ] <- NA
  } else {
    outliers  <- NULL
  }

  if( any( colMeans( is.na(Y) ) == 1 ) )
    stop( 'Entirely blank phenotypes not allowed' )

  if( scale ){
    phenmeans <- apply( Y, 2, mean, na.rm=TRUE )
    phensds   <- apply( Y, 2, sd,   na.rm=TRUE )
    Y0        <- Y
    if( intercept  ){
      Y[,-1]    <- ( ( Y - rep(1,nrow(Y)) %o% phenmeans ) / ( rep(1,nrow(Y)) %o% phensds ) )[,-1]
    } else {
      Y         <- ( Y - rep(1,nrow(Y)) %o% phenmeans ) / ( rep(1,nrow(Y)) %o% phensds )
    }
  }

  if( quantnorm )
    Y <- quantnorm(Y)

  if (missing(K)){
    if( missing( Q ) | missing( lam_K ) )
      stop( 'Either K or its eigen-components, Q and lam_K, must be given' )
    K <- Q %*% diag(lam_K) %*% t(Q)
  }
  if (nrow(K) != nrow(Y)) 
    stop("K and Y are non-conformable")
  if( missing(Q) | missing( lam_K ) ){
    eig.K <- eigen( K, symmetric=TRUE )
    Q     <- eig.K$vec
    lam_K <- eig.K$val
    rm( eig.K )
  }

  #### drop wholly-missing samples
  obs   <- which( rowSums( ! is.na(Y) ) > 0 )
  if( length(obs) == nrow(Y) ){
    eig.K <- list( vec=Q, val=lam_K )
  } else {
    eig.K <- eigen( K[obs,obs], symmetric=TRUE )
  }

  #### run phenix
  vb_out  <- phenix_vb( Y=Y[ obs, ], Q=eig.K$vec, lam_K=eig.K$val, ... )

  #### impute wholly-missing samples
  Yhat            <- Y + NA
  Uhat            <- Y + NA
  Yhat[ obs, ]    <- vb_out$Y
  Uhat[ obs, ]    <- vb_out$U
  if( length( obs ) != nrow(Y) ){
    Yhat[ -obs, ]   <- K[ -obs, obs ] %*% eig.K$vec %*% diag(1/eig.K$val) %*% t(eig.K$vec) %*% vb_out$U
    Uhat[ -obs, ]   <- Yhat[ -obs, ]
  }

  if( scale ){
    if( intercept ){
      Yhat[,-1] <- (  Yhat * ( rep(1,nrow(Y)) %o% phensds ) + rep(1,nrow(Y)) %o% phenmeans )[,-1]
      Uhat[,-1] <- (  Uhat * ( rep(1,nrow(Y)) %o% phensds ) + rep(1,nrow(Y)) %o% phenmeans )[,-1]
    } else {
      Yhat      <-    Yhat * ( rep(1,nrow(Y)) %o% phensds ) + rep(1,nrow(Y)) %o% phenmeans
      Uhat      <-    Uhat * ( rep(1,nrow(Y)) %o% phensds ) + rep(1,nrow(Y)) %o% phenmeans
    }
  }

  if( test ){

    test_out  <- list(
      time      = proc.time()[3] - t0,
      cor_glob  = cor( Yhat[mask], Y.input[mask] ),
      cors      = sapply( test.cols, function(p){
        mask.p  <- intersect( mask, (p-1)*nrow(Y) + 1:nrow(Y) )
        if( length( mask.p ) == 0 ) return( NA )
        cor( Yhat[mask.p], Y.input[mask.p] )
      }),
      mse_glob  = mean( ( Yhat[mask] - Y.input[mask] )^2 ),
      mses      = sapply( test.cols, function(p){
        mask.p  <- intersect( mask, (p-1)*nrow(Y) + 1:nrow(Y) )
        if( length( mask.p ) == 0 ) return( NA )
        mean( ( Yhat[mask.p] - Y.input[mask.p] )^2 )
      })
    )

    return( list( imp=Yhat, U=Uhat, S=vb_out$S, beta=vb_out$beta, h2=vb_out$h2, Q=Q, lam_K=lam_K, vb_out=vb_out, B=vb_out$B, E=vb_out$E, outliers=outliers, test=test_out ) )

  } else {                                                                                                 
    return( list( imp=Yhat, U=Uhat, S=vb_out$S, beta=vb_out$beta, h2=vb_out$h2, Q=Q, lam_K=lam_K, vb_out=vb_out, B=vb_out$B, E=vb_out$E, outliers=outliers ) )
  }

}
