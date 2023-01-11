\name{sim_data}
\alias{sim_data}
\alias{sib_K}
\alias{rpheno}
\alias{rmask}
\title{Simulate phenix data}
\description{Simulates kinship matrix and partially observed phenotype matrix}
\usage{
sib_K( N, fam_size=5 )
rpheno( N, P, K, h2=rep( .5, P ), B, E )
rmask( Y0, miss=.05 )
}
\arguments{
  \item{N}{Number of desired samples}
  \item{fam_size}{Size of each family in kinship matrix}
  \item{P}{Number of desired traits}
  \item{K}{A positive semidefinite kinship matrix describing the relatedness between samples}
  \item{h2}{A length 'P' vector of phenotype heritabilities}
  \item{B}{An optional positive semidefinite matrix of genetic correlations; defaults to Wishart( P, 1/P*diag(P) )}
  \item{E}{An optional positive semidefinite matrix of environmental correlations; defaults to Wishart( P, 1/P*diag(P) )}
  \item{Y0}{A fully-observed matrix of phenotypes}
  \item{miss}{The fraction of entries in Y0 to be set to missing}
}
\value{
  All functions simulate input data for phenix. 'sib_K' returns a block diagonal NxN kinship matrix with blocks corresponding to independent sets of siblings; 'rpheno' returns an NxP matrix of phenotypes; and 'rmask' hides entries in a matrix completely-at-random.
}
\details{ These functions simulate the necessary input data for phenix from the MPMM model; see the MPMM manual page for a definition of this model and the phenix manual page for examples. sib_K generates the matrix of relatedness between samples; rpheno generates heritable traits based on K and the between-trait covariances B and E; rmask then sets some of these entries to missing completely-at-random. }
\author{Andy Dahl}
\examples{
  N       <- 6
  P       <- 2
  (K      <- sib_K(N,fam_size=3))
  (Y      <- rpheno(N,P,K))
  (Ymiss  <- rmask(Y, miss= 0.2))
}

\seealso{\code{\link{MPMM}},\code{\link{phenix}}}
