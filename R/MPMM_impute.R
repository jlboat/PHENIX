MPMM_impute <- function( Y, K, eig.tol=1e-6, standardize=TRUE, ... ){

  if( standardize )
    if( any( abs(apply(Y,2,mean,na.rm=TRUE)) > 1e-8 ) | any( abs( apply(Y,2,sd,na.rm=TRUE) - 1 ) > 1e-8 ) ){
      warning( 'Y unscaled; scaling internally but then returning imputations on original scale' )
      phenmeans <- apply( Y, 2, mean, na.rm=TRUE )
      phensds   <- apply( Y, 2, sd,   na.rm=TRUE )
      Y         <- ( Y - rep(1,nrow(Y)) %o% phenmeans ) / ( rep(1,nrow(Y)) %o% phensds )
    } else {
      standardize <- FALSE
    }

  ##### Eigendecompose (Y-observed parts of) K
  P           <- ncol(Y)
  full_rows   <- which( rowMeans(is.na(Y)) == 0 )
  if( length( full_rows ) < P ){
    warning( paste0( 'Not enough full rows: need', P, ', have ', length( full_rows ), '; returning matrix with 0s replacing NA' ) )
    Y[ is.na(Y) ] <- 0
    return( list( Y=Y ) )
  }

  K.eig       <- eigen( K[ full_rows, full_rows ], symmetric=TRUE )
  if( any( K.eig$val < eig.tol ) ){
    warning(paste( 'Some K eigenvalues are too small; being set to', eig.tol ))
    K.eig$val[ K.eig$val < eig.tol ]  <- eig.tol
  }
  R.val       <- 1/K.eig$val
    
  Y.prime     <- t( K.eig$vec ) %*% Y[ full_rows, ]
  K.svd       <- list( u=K.eig$vec, d=K.eig$val )
  rm( K.eig )

  ####### estimate VCs on full data
  out <- MPMM( Y.prime=Y.prime, Lam.R=R.val, ... )

  ####### impute w/ estimated VCs
  if( ! all( is.na( Y[ - full_rows, ] ) ) ){ ### sporadic missingness
    warning( 'Imputing with MPMM is not advised for sporadically missing phenotypes!' )
    Sigma     <- solve( out$C ) %x% K + solve( out$D ) %x% diag( nrow(Y) )
    o         <- which( ! is.na(c(Y)) )
    X         <- solve( Sigma[ o, o ] ) %*% Y[ o ] ### O( nobs^3 ) = O( N^3 P^3 )
    U         <- matrix( ( solve( out$C ) %x% K )[ ,o ] %*% X, nrow(Y), P )
    Y[-o]     <- Sigma[ -o, o ] %*% X
  } else {
    inv.C.5 <- mat.sqrt( out$C, inv=TRUE )
    U	  <- K[,full_rows] %*% (
      KS_dot_vec( B=Y[full_rows,] %*% solve(inv.C.5) , A=K[full_rows,full_rows], svd.A=K.svd, C=solve(inv.C.5 %*% out$D %*% inv.C.5), inv=TRUE )
    ) %*% inv.C.5
    Y[ -full_rows, ] <- U[ -full_rows, ]
  }

  if( standardize ){
    Y <- Y * ( rep(1,nrow(Y)) %o% phensds ) + rep(1,nrow(Y)) %o% phenmeans
    U <- U * ( rep(1,nrow(Y)) %o% phensds ) + rep(1,nrow(Y)) %o% phenmeans
  }

  return( list( Y=Y, U=U, C=out$C, D=out$D, VC_out=out ) )

}

MPMM <- function( Y, K, Y.prime, Lam.R,
          tol=1e-4, maxit=1e4, warm.updates=TRUE, trace=FALSE, print.like=trace,
          lmm_init=TRUE, start.C, start.D, eig.tol=1e-6, standardize=missing(Y.prime)
){

  if( missing( Y.prime ) | missing( Lam.R ) ){
    if( any( is.na( Y ) ) )
      stop( 'Missing phentoypes not allowed' )
    if( standardize )
      if( any( abs(apply(Y,2,mean,na.rm=TRUE)) > 1e-8 ) | any( abs( apply(Y,2,sd,na.rm=TRUE) - 1 ) > 1e-8 ) ){
        warning( 'Y unscaled; scaling now' )
        phenmeans <- apply( Y, 2, mean, na.rm=TRUE )
        phensds   <- apply( Y, 2, sd,   na.rm=TRUE )
        Y         <- ( Y - rep(1,nrow(Y)) %o% phenmeans ) / ( rep(1,nrow(Y)) %o% phensds )
      } else {
        standardize <- FALSE
      }
    if( trace )
      print( 'eigendecomposing K...' )
    K.eig       <- eigen( K, symmetric=TRUE )
    if( any( K.eig$val < eig.tol ) ){
      warning(paste( 'Some K eigenvalues are too small; being set to', eig.tol ))
      K.eig$val[ K.eig$val < eig.tol ]  <- eig.tol
    }
    Lam.R       <- 1/K.eig$val
    Y.prime     <- t( K.eig$vec ) %*% Y
    rm( K.eig, Y, K )
    if( trace )
      print('preprocessing done')
  } else {
    if( standardize )
      warning( 'Skipping standardization; using Yprime as-is instead' )
  }

  miss.VC <- ( missing( start.C ) | missing( start.D ) )
  null.VC <- ifelse( miss.VC, TRUE, is.null( start.C ) | is.null( start.D ) )

  if( miss.VC | null.VC ){
    if( lmm_init ){
      lmms  <- lapply( 1:ncol(Y.prime), function(p) fit_lmm(yprime=Y.prime[,p],Lam.K=1/Lam.R) )
      C     <- diag( sapply( lmms, function( lmm ) 1/lmm$sig_g ) )
      D     <- diag( sapply( lmms, function( lmm ) 1/lmm$sig_e ) )
    } else {
      C   <- 1/2*diag(ncol(Y.prime))
      D   <- 1/2*diag(ncol(Y.prime))
    }
  } else {
    C     <- start.C
    D     <- start.D
  }
  C.out <- list( par=C, invpar=solve(C) )
  D.out <- list( par=D, invpar=solve(D) )

  del  	<- numeric()
  times <- numeric(7)
  if(trace){
    cat( "\n", sprintf( "%4s, %8s, %8s, %8s, %10s, %10s, %10s", "iter", "delta", "delta C", "delta D", "Om time", "C  time", "D  time" ) )
    if( print.like )
      cat( sprintf( ", %10s", "-loglike" ) )
    cat("\n")

    ###### initialization likelihood
    if( print.like ){
      like <- MPMM_like( Y.prime=Y.prime, Lam.R=Lam.R, C=C, D=D )
      ll_path <- numeric()
      cat( sprintf( "%5d, %8s, %8s, %8s, %10s, %10s, %10s, %10.1f", 0, NA, NA, NA, 0, 0, 0, -like ), "\n" )
    } else {
      ll_path <- NA
    }
  } else {
    ll_path <- NA
  }

  for( t in 1:maxit ){
    
    oldD <- D
    oldC <- C

    times[1]  <- times[1] + system.time(
      Omegas  <- MPMM_Estep( C, D, Lam.R, Y.prime )
    )[3]
      
    times[2]  <- times[2] + system.time({
      C.out   <- MPMM_Mstep( Omega=Omegas$C, pars=list(type=NA,rho=NA), warm=warm.updates, init=C.out )
      C       <- C.out$par
    })[3]

    times[3]  <- times[3] + system.time({
      D.out   <- MPMM_Mstep( Omega=Omegas$D, pars=list(type=NA,rho=NA), warm=warm.updates, init=D.out )
      D       <- D.out$par
    })[3]

    rmse      <- function(x,y) sqrt( mean( (x-y)^2 ) )
    del.C     <- rmse(C,oldC) / rmse(C[upper.tri(C)],0+1e-10) ### offset for numerical stability
    del.D     <- rmse(D,oldD) / rmse(D[upper.tri(D)],0+1e-10) ### offset for numerical stability
    del[t]    <- max( del.C, del.D )

    if(trace){
      cat( sprintf( "%5d, %8.4f, %8.4f, %8.4f, %10.4f, %10.4f, %10.4f", t, del[t], del.C, del.D, times[1], times[2], times[3] ) )
      if( print.like ){
        old_like  <- like
        ll_path[t]<- like
        like      <- MPMM_like( Y.prime=Y.prime, Lam.R=Lam.R, C=C, D=D )
        cat( sprintf( ", %9.1f", -like ) )
        #if( like + 1e-8 < old_like )
        #  warning(paste('Likelihood decreased by:', like - old_like ))

        if( (like-old_like)/like < -1e-6 )
          warning(paste('Likelihood decreased by:', like - old_like ))
      }
      cat("\n")
    }
    if( del[t] < tol ) break  
  }

  if( standardize ){
    C.unscaled  <- diag( 1/phensds ) %*% C %*% diag( 1/phensds )
    D.unscaled  <- diag( 1/phensds ) %*% D %*% diag( 1/phensds )
  } else {
    C.unscaled  <- NA
    D.unscaled  <- NA
  }
  like      <- MPMM_like( Y.prime=Y.prime, Lam.R=Lam.R, C=C, D=D )
  h2        <- diag(solve(C)) / ( diag(solve(C))+diag(solve(D)) )
  return( list( del=del, h2=h2, C=C, D=D, times=times, niter=t, conv=( t!=maxit ), C.unscaled=C.unscaled, D.unscaled=D.unscaled, ll=like, ll_path=ll_path ) )
}

MPMM_Estep <- function( C, D, Lam.R, Y.prime ){
  
  N	<- nrow(Y.prime)
  P	<- ncol(Y.prime)

  inv.sqrt.C <- mat.sqrt(C,inv=TRUE)
  inv.sqrt.D <- mat.sqrt(D,inv=TRUE)    

  svd.1 <- svd( inv.sqrt.C %*% D %*% inv.sqrt.C )
  svd.2 <- svd( inv.sqrt.D %*% C %*% inv.sqrt.D )

  Sig.1	<- inv.sqrt.C %*%  trp_KS( svd.keep = svd.1, lam.kill=Lam.R ) %*% inv.sqrt.C
  Sig.2	<- inv.sqrt.D %*%  trp_KS( svd.keep = svd.2, lam.kill=1/( Lam.R ) ) %*% inv.sqrt.D

  M	<- KS_dot_vec( B=Y.prime %*% D %*% inv.sqrt.C, svd.A=list( d=Lam.R ), svd.C=svd.1, inv=TRUE, diag.A=TRUE ) %*% inv.sqrt.C

  M.1	<- t(Y.prime-M) %*% (Y.prime-M)
  M.2	<- t(M) %*% ( matrix(Lam.R,N,P) * M )
	
  return( list( D = 1/N * ( Sig.1 + M.1 ), C = 1/N * ( Sig.2 + M.2 ) ) )
  
}

MPMM_Mstep  <- function( Omega, pars, warm, init ){
  P <- dim(Omega)[1]
  if( is.na( pars$type ) ){
    par			<- solve( Omega )
  }
  return( list( par = par ) )  
}

MPMM_like <- function( Y.prime, Lam.R, C, D ){
    
  C.5     <- mat.sqrt(C)
  eig.Delta <- eigen( C.5 %*% solve(D) %*% C.5,symmetric=TRUE )

  Lambda  <- rep(1,ncol(Y.prime)) %x% (1/Lam.R)+ eig.Delta$val %x% rep(1,nrow(Y.prime))
  x       <- c( Y.prime %*% C.5 %*% eig.Delta$vec )

  lik2    <- nrow(Y.prime)*logdet(C) - sum( log( Lambda ) ) - sum(( vec.pinv( sqrt( Lambda ) ) * x )^2)
  return( lik2/2 )
    
}
