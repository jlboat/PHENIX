\name{LMM}

\alias{LMM}
\alias{fit_lmm}
\alias{lmm_obj}
\alias{lmm.impute}

\title{Fit a random effect model}
\description{Fit univarite random effect models and predict missing phenotypes and breeding values.}

\usage{
  lmm.impute( Y, K, eig.marg, trace=FALSE )
  fit_lmm( y, yprime, K, eig.K, Lam.K, constrained=TRUE )
  lmm_obj( delta, yprime, Lam.K )
}
\arguments{
  \item{Y}{Centered and scaled matrix of phenotypes}
  \item{K}{ A positive semidefinite kinship matrix describing the genetic covariance between individuals }
  \item{eig.marg}{ An optional list containing eigendecompositions of K; entry p eigendecomposes K for the individuals on which phenotype p is observed }
  \item{trace}{Prints algorithm's progress through loop over traits if true}
  \item{y}{Centered and scaled vector of phenotypes}
  \item{yprime}{Kinship-whitened vector of phenotypes}
  \item{eig.K}{ The eigendecomposition of K for individuals with observed phenotypes }
  \item{Lam.K}{Eigenvalues of K for individuals with observed phenotypes }
  \item{delta}{Ratio of environmental and genetic variance components}
  \item{constrained}{Whether fitted heritability is constrained to [0,1]}
}
\value{

  'lmm.impute' returns imputed phenotypes, predicted breeding values and variance components for all phenotypes:
    \item{Yhat}{Imputed phenotype matrix}
    \item{h2}{Vector of heritibailites for each phenotype}
    \item{U}{Predicted matrix of breeding values}

  'fit_lmm' returns variance components for a single phenotype:
    \item{sig_g}{Genetic variance component}
    \item{sig_e}{Environemntal variance component}
    \item{h2}{Phenotype heritability}

  'lmm_obj' returns the log-likelihood.
}

\author{Andy Dahl}


\details{ 'fit_lmm' fits sig_g and sig_e by maximum likelihood using the model

y ~ N( 0, sig_g^2 K ) + N( 0, sig_e^2 I )

where y is a fully observed vector of phenotypes and K is a known kinship matrix. yprime is a version of y that has been whitened by the eigenvectors of K (i.e. if K = U diag(1/Lam.R) U^T is an eigendecomposition, yprime = U^T y) so that entries of yprime are independent.

'lmm.impute' uses the variance components from univariate lmms (on the observed phenotypes) to predict missing entries of Y and a matrix of breeding values U via BLUPs (best linear unbiased predictions).
}

\examples{

  # simulate data
  N   <- 500
  P   <- 10
  K   <- sib_K( N=N )             # indepedendent sets of sibs
  Y0  <- rpheno( N=N, P=P, K=K )  # simulate phenotype matrix
  Y   <- rmask( Y0, miss=.05 )    # hide phenotypes

  # run an lmm on one trait
  fit_lmm(y=Y[,1], K=K)

  # create+decompose K for observed samples at each phenotype
  eig.marg  <- lapply( 1:ncol(Y), function(p){
    obs   <- which( !is.na(Y[,p]) )
    K_loc <- K[ obs, obs ]
    return( eigen( K_loc ) )
  })
    
  # fit all lmms and predict phenotypes + breeding values
  out   <- lmm.impute( Y, K, eig.marg=eig.marg, trace=TRUE )
  out$h2            # estimated heritabilies
  head( out$Yhat )  # imputed phenotypes
  
}
\keyword{ algebra }
