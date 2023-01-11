sib_K <- function( N, fam_size=5 ){

  if( N != abs(round( N )) )
    stop( 'N must be a positive integer' )

  if( fam_size != abs(round( fam_size )) )
    stop( 'fam_size must be a positive integer' )

  num_fams  <- N / fam_size
  if( num_fams != round( num_fams ) )
    stop( 'fam_size must divide N' )

  K   <- matrix( 0, N, N )
  for( family in 1:num_fams ){
    fam_inds    <- (family-1)*fam_size + 1:fam_size
    K[ fam_inds, fam_inds ] <- .5
  }
  diag(K) <- 1
  K
}

rpheno  <- function( N, P, K, h2=rep( .5, P ), B, E ){

  if( ! all( dim(K) == N ) )
    stop( 'K must be a NxN matrix' )

  if( missing(B) ){
    B <- rWishart( 1, P, 1/P*diag(P) )[,,1]
  } else {
    if( ! all( dim(B) == P ) )
      stop( 'B must be a PxP matrix' )
  }

  if( missing(E) ){
    E <- rWishart( 1, P, 1/P*diag(P) )[,,1]
  } else {
    if( ! all( dim(E) == P ) )
      stop( 'E must be a PxP matrix' )
  }

  B <- diag(sqrt(h2)  ) %*% cov2cor(B) %*% diag(sqrt(h2)  )
  E <- diag(sqrt(1-h2)) %*% cov2cor(E) %*% diag(sqrt(1-h2))

  Y <- mat.sqrt(K) %*% matrix( rnorm(N*P), N, P ) %*% mat.sqrt(B) +
       matrix( rnorm(N*P), N, P ) %*% mat.sqrt(E)
  Y

}

rmask <- function( Y0, miss=.05 ){

  if( miss < 0 | miss > 1 )
    stop( 'miss must be between 0 and 1' )

  mask  <- sample( length(Y0), round( length(Y0) * miss ) )
  Y     <- Y0
  Y[mask] <- NA
  Y

}
