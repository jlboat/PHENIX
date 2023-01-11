mv_gwas <- function( Y, G, K, Y.prime, G.prime, lam_K, trace=FALSE, eig.tol=1e-6, ...
){

  if( missing( Y.prime ) | missing( G.prime ) | missing( lam_K ) ){

    if( any( is.na( Y ) ) )
      stop( 'Partially-missing Y not supported' )

    if( trace )
      print( 'Eigendecomposing K...' )
    K.eig       <- eigen( K, symmetric=TRUE )
    if( any( K.eig$val < eig.tol ) ){
      warning(paste( 'Some K eigenvalues are too small; being set to', eig.tol ))
      K.eig$val[ K.eig$val < eig.tol ]  <- eig.tol
    }
    lam_K       <- K.eig$val

    if( trace )
      print( 'Whitening Y...' )
    Y.prime <- t( K.eig$vec ) %*% Y             ### whiten phenos 

    if( trace )
      print( 'Whitening G...' )
    G.prime <- t( K.eig$vec ) %*% scale( G )    ### whiten genos 

    rm( K.eig, Y, K, G )
  }
  
  ### fit C, D under null
  if( trace )
    print( 'Running null....' )
  null  <- MPMM( Y.prime=Y.prime, Lam.R=1/lam_K, trace=trace, ... )

  ### compute alt lls using C, D from null
  if( trace )
    print( 'Running LRTs...' )
  alt_ll  <- alt_ll_fxn( Y.prime, G.prime, null$C, null$D, lam_K=lam_K, trace=trace )

  ### compute LRT
  lrt     <- -2*( null$ll - alt_ll )
  p       <- -pchisq( lrt, ncol(Y.prime), log.p=TRUE, lower.tail=FALSE ) / log( 10 )

  return( list( lrt=lrt, p=p, null=null ) )

}

alt_ll_fxn  <- function( Y.prime, G.prime, C, D, lam_K, trace=FALSE
){

  N   <- nrow(Y.prime)
  P   <- ncol(Y.prime)

  #### lots of O(P^3)
  C.5   <- mat.sqrt( C )
  eig.P <- eigen( C.5 %*% solve(D) %*% C.5 )
  lam.P <- eig.P$val
  L     <- C.5 %*% eig.P$vec
  Linv  <- solve( L )
  rm( eig.P )

  #### O(N P^2)
  Y.trans <- ( Y.prime %*% L )

  #### O( NP + P^3 )
  Lambda  <- rep(1,P) %x% lam_K + lam.P %x% rep(1,N)
  Z       <- Y.trans / matrix( Lambda, N, P )

  l_base  <- N*logdet(C) - sum(log( Lambda ))
  sqrt.Lambda <- sqrt( vec.pinv( Lambda ) )

  #### O( ncol(G.prime) [ NP + P^2 ] )
  sapply( 1:ncol(G.prime), function(l){
    g <- matrix(G.prime[,l],N,1)
    if( trace )
      cat( 'Locus #: ', l, '\n' )

    ### uses C, D from null to maximize alt
    omega   <- 1/sapply( 1:P, function(p) t(g) %*% diag( 1/Lambda[ (p-1)*N + 1:N ] ) %*% g )
    betahat <- ( t(g) %*% Z ) %*% ( ( omega %o% rep( 1, P ) ) * Linv )

    ### compute alt ll using C, D from null
    x       <- c( Y.trans - g %*% ( betahat %*% L ) )
    return( ll=1/2*( l_base - sum( (sqrt.Lambda*x)^2 ) ) )
  })

}
