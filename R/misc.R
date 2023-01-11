mat.sqrt <- function(X,inv=FALSE){
  if( length(X) == 1 )
    if( inv ){
      return(1/sqrt(X))
    } else {
      return(sqrt(X))
    }
  eig <- eigen(X,symmetric=TRUE)
  Q   <- eig$vec

  negs  <- which( eig$val < 0 )
  if( length( negs ) > 0 ){
    warning( paste( 'There are', length( negs ), 'negative eigenvalues; they are being set to 1e-8' ) )
    eig$val[ negs ] <- 1e-8
  }
  
  if( ! inv ){
    Lam <- diag( sqrt(eig$val) )
  } else {
    Lam <- diag( 1/sqrt(eig$val) )
  }
  
  return( Q %*% Lam %*% t(Q) )
}

schur  <- function( X, m ){

  if( min(m) < 0 | max(m) > ncol(X) )
    stop( 'm must have entries in 1, 2, ..., ncol(X)' )

  X[m,m] - X[m,-m,drop=FALSE]%*%solve(X[-m,-m])%*%t(X[m,-m,drop=FALSE])
}

vec.pinv  <- function( x, tol=1e-8 )
  as.numeric( sapply( x, function(y) ifelse( abs(y) > tol, 1/y, 0 ) ) )

quantnorm <- function(Y){

  P <- ncol(Y)
  for( p in 1:P ){

    if( all( Y[,p] == 1 ) )
      next

    obs       <- which( ! is.na(Y[,p]) )
    ranks     <- rank( Y[obs,p] )
    quantiles <- qnorm( ranks/(length(ranks)+1) )
    Y[obs,p]  <- quantiles

  }
  Y   <- scale( as.matrix(Y) )
  Y   <- matrix( Y, nrow(Y), ncol(Y) )
  return( Y )
}

logdet <- function(X)
  ifelse(length(X) == 1, log(X), as.numeric(determinant(X, logarithm = TRUE)$mod))
