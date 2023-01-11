phenix_vb <- function( Y, Q, lam_K, N=nrow(Y), P=ncol(Y), M=min(c(N,P)),
    tau=0, e=P+5, E.inv=solve( diag(P)/(e-P-1) ),
    trace=0, reltol=1e-8, maxit=1e3, cutoff=1e-3, ... )
{

###### check eigenvalues of K
  if( any( lam_K < 0 ) ){
    warning ( "K has negative eigenvalues:")
    print   ( c( "range:", range( lam_K[ lam_K < 0 ] ) ) )
    print   ( c( "number:", length( lam_K[ lam_K < 0 ] ) ) )
    print   ( 'Setting negative eigenvalues to 0' )
    lam_K[ lam_K < 0 ]  <- 0
  }

###### Set hyperparameters and constants
  N 	  <- nrow(Y)
  P 	  <- ncol(Y)
  eprime  <- e + N  
 
  fullInds	<- which( rowMeans(is.na(Y)) == 0 )  
  Y.indic	<- is.na(Y) + 0
  j		<- 0
  done		<- rep( 0, N )
  if( length( fullInds ) > 0 )
    done[fullInds]	<- 1
    
  j2ns	<- j2miss	<- list()
  j2len	<- numeric()
  while( mean(done) < 1 ){
    j			<- j+1
    n			<- which( done == 0 )[1]
    miss.pattern	<- Y.indic[n,]  	
    row.indic		<- apply( Y.indic, 1, function(y) mean( ( y + miss.pattern ) %% 2 ) )
    chosen		<- which( row.indic == 0 )
    j2ns[[ j ]]		<- chosen
    j2len[ j ]		<- length(chosen)
    j2miss[[ j ]]	<- which( miss.pattern == 1 )		
    done[chosen ]	<- 1
  }
  n.miss.types	<- j

  if( trace > 0 )
    print( "starting initialization...")

####### Initialize Y and Omega
  if( all( Y[,1] == 1 ) ){
    if( trace )
      print('intercept found')
    out      	<- MVN_impute( Y[ ,-1 ], intercept=FALSE, ... )
    out$Y       <- cbind( 1, out$Y )
    out$Sigma   <- cbind( c( 1, rep(0,P-1) ), rbind( 0, out$Sigma ) )
  } else {
    out      	<- MVN_impute( Y, intercept=FALSE, ... )
  }
  mu_Y    	<- out$Y
  Omega.inv	<- out$Sigma
  Omega		<- solve(Omega.inv)
  svd.Omega	<- svd(Omega,LINPACK=TRUE)
 
  if( n.miss.types == 0 ){
    Sigma_Y	<- array( 0, dim=c( 2, P, P ) )
  } else {
    Sigma_Y	<- array( 0, dim=c(n.miss.types,P,P) )
    for( j in 1:n.miss.types ){
      m	    	        <- j2miss[[j]]
      Sigma_Y[j,m,m]	<- j2len[j]	* solve( Omega[m,m] )    
    }
  }
  
####### Initialize S and beta
  out		<- svd( mu_Y, LINPACK=TRUE )
  mu_S		<- ( out$u %*% diag(sqrt(out$d))    )[ , 1:M ]
  mu_beta	<- ( diag(sqrt(out$d)) %*% t(out$v) )[ 1:M,  ]
   
  V_beta	  <- V_S            <- diag(M)
  svd.V_beta	  <- svd.V_S	  <- list( u=diag(M), d=rep(1,M) )
  svd.Omega_beta  <- svd.Omega
  	
################ Set up looping
  counter	      <- 1
  lik.path	      <- zero.path    <- numeric()
  lik.path[counter]   <- lik	      <- phenix_likelihood( Y, Q, lam_K, N, P, M, tau, E.inv,
    mu_Y, Sigma_Y, Omega, mu_beta, mu_S, eprime,
    svd.Omega, svd.V_S, svd.V_beta, svd.Omega_beta,
    j2miss, j2len, n.miss.types )
  zero.path[counter]  <- zeros	      <- sum( sapply( 1:M, function(m) mean(abs(mu_S[,m]%o%mu_beta[m,])) ) < cutoff )
  if( trace > 0 ){
    print( "initialization done")
    cat( sprintf("%10s, %10s, %10s, %10s", "counter", "loglik", "Zeroes", "delta l"), "\n" )  
    cat( sprintf("%10d, %10.3e, %10d, %10.3e \n", counter, lik, zeros, NA) )
  }
  repeat{

    ############### beta
    Omega_beta 	    <- Omega
    svd.Omega_beta  <- svd.Omega
    V_S		    <- t(mu_S) %*% mu_S + trp_KS( svd.keep=svd.V_beta, lam.kill=1/lam_K )
    svd.V_S	    <- svd( V_S, LINPACK=TRUE )
    mu_beta	    <- KS_dot_vec( svd.C = list( u=svd.Omega$u, d=tau/svd.Omega$d ), svd.A = svd.V_S, B=t(mu_S) %*% mu_Y, inv=TRUE )
    ############### S
    V_beta	    <- mu_beta %*% Omega %*% t(mu_beta) + trp_KS( lam.kill = tau/svd.Omega$d, svd.keep=svd.V_S )
    svd.V_beta	    <- svd(V_beta,LINPACK=TRUE)
    mu_S	    <- KS_dot_vec( svd.C = svd.V_beta, svd.A = list( u = Q, d = 1/lam_K ), B = mu_Y %*% Omega %*% t(mu_beta), inv=TRUE )

    ############### beta, again
    V_S		    <- t(mu_S) %*% mu_S + trp_KS( svd.keep=svd.V_beta, lam.kill=1/lam_K )
    svd.V_S	    <- svd( V_S )
    mu_beta	    <- KS_dot_vec( svd.C = list( u=svd.Omega$u, d=tau/svd.Omega$d ), svd.A = svd.V_S, B=t(mu_S) %*% mu_Y, inv=TRUE )

    ###################################### Update Omega  	
    mu_Y.d	    <- mu_Y - mu_S %*% mu_beta	
    trp.s	    <- trp_KS( svd.keep = svd.V_beta, lam.kill = 1/lam_K )
    if( tau != 0 ){
      Delta	    <- trp_KS( svd.keep=svd.Omega_beta, lam.kill=tau/svd.V_S$d ) 			
    } else {
      Delta	    <- solve( Omega_beta ) * M
    }
    Omega.inv	    <- 1/eprime * ( t(mu_Y.d)%*%mu_Y.d + apply( Sigma_Y, 2:3, sum ) + t(mu_beta)%*%trp.s%*%mu_beta + Delta + E.inv )
    svd.Omega.inv   <- svd(Omega.inv,LINPACK=TRUE)
    svd.Omega	    <- list( u = svd.Omega.inv$u, d=1/svd.Omega.inv$d )
    Omega	    <- svd.Omega$u %*% diag(svd.Omega$d) %*% t(svd.Omega$u)
    ################################### Update Y  	
    U       <- mu_S %*% mu_beta
    mu_Y    <- Y
    if( n.miss.types != 0 ){
      Sigma_Y		<- array( 0, dim=c(n.miss.types,P,P) )
      for( j in 1:n.miss.types ){
        m	    	<- j2miss[[j]]
        o  		<- (1:P)[-m]	    
        ns		<- j2ns[[j]]
        mu_Y[ns,m]	<- U[ns,m] + ( Y - U )[ns,o,drop=FALSE]%*%solve(Omega.inv[o,o])%*%Omega.inv[o,m,drop=FALSE]
        Sigma_Y[j,m,m]	<- j2len[j]	* solve( Omega[m,m] )    
      }
    }
	
    #### compute lik, zeros, update counter
    counter             <- counter+1   
    oldlik		<- lik  
    lik.path[counter]   <- lik		<- phenix_likelihood( Y, Q, lam_K, N, P, M, tau, E.inv,
    mu_Y, Sigma_Y, Omega, mu_beta, mu_S, eprime,
    svd.Omega, svd.V_S, svd.V_beta, svd.Omega_beta,
    j2miss, j2len, n.miss.types )
    zero.path[counter]  <- zeros	<- sum( sapply( 1:M, function(m) mean(abs(mu_S[,m]%o%mu_beta[m,])) ) < cutoff )  
    delta 		<- (lik - oldlik)/abs(lik)

    if( is.na(lik) ){
      save.image("/data/emu/not-backed-up/dahl/na_lik.Rdata")
      stop("na likelihood")
    }
    #### print, determine whether to break
    if( trace > 0 )
      cat( sprintf("%10d, %10.3e, %10d, %10.3e \n", counter, lik, zeros, delta) )
    if( delta < -reltol )
      if( abs( delta ) < 1e-12 ){ ### neighborhood of numerical error
        warning(paste("KL decreased: delta ll =", delta ))
      } else {
        stop(   paste("KL decreased: delta ll =", delta ))
      }
    if ( delta < reltol ){
      break
    } else if ( counter > maxit ){
      warning(paste("Failed Convergence: delta =",delta))
      break
    }

  }

  B   <- t(mu_beta) %*% mu_beta
  E   <- solve( (eprime-P-1)/eprime * Omega )
  h2  <- diag(B) / ( diag(B) + diag(E))

  vb_pars	<- list( S=mu_S, beta=mu_beta, Omega=Omega, Sigma_Y=Sigma_Y, V_beta=V_beta, V_S=V_S, eprime=eprime, tau=tau, M=M, Q=Q, E.inv=E.inv  )

  return( list(
    Y         = mu_Y, 
    U         = U,
    S         = mu_S,
    beta      = mu_beta,
    counter   = counter, 
    lik.path  = lik.path, 
    zero.path = zero.path, 
    lik       = lik,
    h2        = h2,
    B         = B,
    E         = E,
    vb_pars   = vb_pars,
    tau       = tau, 
    fitted.M  = P - zeros
  ) )

}

phenix_likelihood <- function( Y, Q, lam_K, N, P, M, tau, E.inv,
    mu_Y, Sigma_Y, Omega, mu_beta, mu_S, eprime,
    svd.Omega, svd.V_S, svd.V_beta, svd.Omega_beta,
    j2miss, j2len, n.miss.types ){
  
  if( n.miss.types == 0 ){
    log.det.Y	<- 0
  } else {
    log.det.Y  = rep( 0, n.miss.types )
    for( j in 1:n.miss.types ){
      miss	      <- j2miss[[j]]
      log.det.Y[ j ]  <- j2len[j] * ( logdet( Sigma_Y[ j, miss, miss ] ) - length(miss) * log(j2len[j]) )
    }
  }
  lam_V_s	<- svd.V_S$d
  lam_V_beta	<- svd.V_beta$d  
  lam_omega	<- svd.Omega_beta$d

  if( any( lam_K == 0 ) )
    lam_K[ lam_K == 0 ] <- 1e-6
  	
  mu_Y.d	<- mu_Y - mu_S %*% mu_beta	

  trp.s		<- trp_KS( svd.keep=svd.V_beta, lam.kill= 1/lam_K )

  Delta		<- trp_KS( svd.keep = svd.Omega_beta, lam.kill = tau/lam_V_s )

  mu_S.term     <- sum( ( matrix( 1/sqrt(lam_K), N, M ) * ( t(Q) %*% mu_S ) )^2 )
		
  ll 		<- 1/2*(
                  sum( log.det.Y ) + eprime * sum( log( svd.Omega$d ) ) - tau*( sum(mu_beta^2) + sum( 1 / ( lam_omega %x% lam_V_s + rep(tau,M*P) ) ) ) - sum( log( lam_omega %x% lam_V_s + rep(tau,M*P) ) ) 
                  - mu_S.term - sum( 1 / ( lam_V_beta %x% lam_K + rep(1,N*M) ) ) - sum( log( lam_V_beta %x% rep(1,N) + rep(1,M) %x% (1/lam_K) ) )  
                  - tr( Omega %*% ( t(mu_Y.d)%*%mu_Y.d + apply( Sigma_Y, 2:3, sum ) + t(mu_beta) %*% trp.s %*% mu_beta + Delta + E.inv ) ) 								
  )

  return( ll )
}
