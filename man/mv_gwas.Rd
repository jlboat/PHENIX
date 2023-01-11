\name{mv_gwas}

\alias{mv_gwas}
\alias{alt_ll_fxn}

\title{Perform multiphenotype GWAS via approximate LRT}

\usage{
mv_gwas( Y, G, K, Y.prime, G.prime, lam_K, trace=FALSE, eig.tol=1e-6, ... )
alt_ll_fxn( Y.prime, G.prime, C, D, lam_K, trace=FALSE )
}

\arguments{
  \item{Y}{Fully observed phenotype matrix}
  \item{G}{Genotype matrix with SNPs as columns}
  \item{K}{A positive semidefinite kinship matrix describing the genetic covariance between individuals}
  \item{Y.prime}{Kinship pre-whitened phenotype matrix}
  \item{G.prime}{Scaled, then kinship pre-whitened genotype matrix}
  \item{lam_K}{Eigenvalues of the kinship matrix K}
  \item{trace}{Prints progress through the algorithm if true}
  \item{eig.tol}{Eigenvalues of K below eig.tol are threshholded to 0}
  \item{...}{Arguments passed to MPMM, used to fit VCs under the null}
  \item{C,D}{Inverse variance components from fitting MPMM under the null}
}
\value{

'alt_ll_fxn' returns the (approximate) maximum likelihood under the alternatives. 

'mv_gwas' returns a list:
  \item{lrt}{Vector of LRT statistics for each SNP in G}
  \item{p}{Vector of LRT p-values for each SNP in G}
  \item{null}{Result of running 'MPMM' under the null}

}

\author{Andy Dahl}

\description{Approximate multivariate GWAS method to test whether each of L SNPs is associated with (at least one of) the P traits.}

\details{ For each vector of SNPs g, or column of G, these functions fit the model

Y ~ g beta + MN( 0, K, C^{-1} ) + MN( 0, I, D^{-1} )

where Y is a fully observed matrix of phenotypes, K is a known kinship matrix, and C and D are precisions describing the heritable and environmental correlations between phenotypes, respectively. The LRTs test the hypothesis that beta != 0; the approximation is that C and D are fit only once, under the null that beta=0, rather than fit jointly with beta at each SNP.

The algorithm requires only Y.prime, G.prime and lam_K; if any of these are missing, however, then Y, G and K must all be provided so the former can be computed. 
}

\examples{

  # dummy data
  N   <- 200
  P   <- 3
  Y   <- matrix( rnorm(N*P), N, P )
  y   <- c(Y)
  g   <- as.numeric(scale(rbinom( N, 2, .5 )))
  K   <- rWishart( 1, N, diag(N) )[,,1]
  C   <- rWishart( 1, P, diag(P) )[,,1]
  D   <- rWishart( 1, P, diag(P) )[,,1]
  eig.K <- eigen(K)
  Q     <- eig.K$vec
  lam_K <- eig.K$val

  # fast log-like
  ll_fast <- alt_ll_fxn( Y.prime=t(Q)\%*\%Y, G.prime=t(Q)\%*\%g, C, D, lam_K=lam_K )

  # check against naive computation: 
  # first compute betahat then evaluate ll

  # two slow ways to compute beta:
  # first:
  Sigma <- solve(C) \%x\% K + solve(D) \%x\% diag(N)
  trans <- mat.sqrt( Sigma, inv=TRUE )
  X     <- diag(P) \%x\% matrix(g,N,1)
  yt    <- trans \%*\% y
  Xt    <- trans \%*\% X
  betahat1 <- as.numeric( coef( lm( yt ~ Xt-1 ) ) )

  # second:
  Omega <- solve(Sigma)
  betahat  <- as.numeric( solve( t(X) \%*\% Omega \%*\% X ) \%*\% t(X) \%*\% Omega \%*\% y )
  all.equal( betahat, betahat1 )

  # ll evaluation using betahat
  ldet  <- -logdet( Sigma )
  resid <- c( Y - g \%o\% betahat )
  quad  <- -as.numeric( t( resid ) \%*\% solve( Sigma ) \%*\% resid )
  all.equal( ll_fast, ll_slow  <- 1/2*( ldet + quad ) )

  # run multiphenotype 'GWAS' on 2 simulated SNPs
  G   <- cbind( g, scale(rbinom( N, 2, .5 )) )  # 'genome-wide'
  out <- mv_gwas( Y=Y, G=G, K=K, maxit=1e3 )    # 'GWAS'
  out$lrt
  out$null$h2

}

\seealso{\code{\link{MPMM}}}

\references{ Andrew Dahl,  Valentina Iotchkova,  Amelie Baud,  Asa Johansson,  Ulf Gyllensten, Nicole Soranzo, Richard Mott, Andreas Kranis \& Jonathan Marchini (2016) A multiple-phenotype imputation method for genetic studies, \url{http://www.nature.com/ng/journal/v48/n4/full/ng.3513.html}, Nature Genetics, 48(4), 466-472 }
