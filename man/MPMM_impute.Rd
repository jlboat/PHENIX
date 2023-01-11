\name{MPMM_impute}
\alias{MPMM_impute}
\title{Phenotype imputation with MultiPhenotype Mixed Models}
\details{Fit an MPMM on fully-observed samples then impute via BLUPs}

\usage{
MPMM_impute( Y, K, eig.tol=1e-6, standardize=TRUE, ... )
}
\arguments{
  \item{Y}{\# samples x \# phenotypes matrix of partially observed phenotypes}
  \item{K}{kinship matrix of relatedness across samples}
  \item{eig.tol}{Eigenvalues of K below eig.tol are thresholded to 0}
  \item{standardize}{If TRUE and Y is unscaled, a warning is produced and Y is scaled internally}
  \item{\dots}{Arguments passed to MPMM, which is used to fit variance components for BLUPs}
}
\value{
  \item{Y}{Imputed phenotype matrix}
  \item{U}{Matrix of estimated breeding values}
  \item{C}{Inverse genetic variance component, estimated by MPMM, used to impute Y}
  \item{D}{Inverse environmental variance component, estimated by MPMM, used to impute Y}
  \item{VC_out}{Output from MPMM run on fully-observed rows of Y}
}

\examples{

  # simulation parameters
  N   <- 300
  P   <- 4
  K   <- sib_K( N=N )             # indepedendent sets of sibs
  Y0  <- rpheno( N=N, P=P, K=K )  # simulate phenotype matrix
  Y   <- rmask( Y0, miss=.05 )    # hide phenotypes

  # run MPMM_impute
  out <- MPMM_impute(Y,K)
  head( out$Y )

}

\description{The method first fits an MPMM on fully-observed samples (rows of Y). Using the resulting variance components, missing entries of Y are imputed to their conditional expectation given the observed entries. If Y is unscaled (has nonzero mean or non-unit variance), it is scaled internally but imputed values are returned on the original scale.}

\author{ Andy Dahl }
\note{The imputation step is O(N^3 P^3) and, when N is large, can be prohibitively computationally expensive even when fitting the MPMM is cheap.}

\references{ Andrew Dahl,  Valentina Iotchkova,  Amelie Baud,  Asa Johansson,  Ulf Gyllensten, Nicole Soranzo, Richard Mott, Andreas Kranis \& Jonathan Marchini (2016) A multiple-phenotype imputation method for genetic studies, \url{http://www.nature.com/ng/journal/v48/n4/full/ng.3513.html}, Nature Genetics, 48(4), 466-472 }

\seealso{\code{\link{MPMM}}, \code{\link{MPMM_helpers}}}
