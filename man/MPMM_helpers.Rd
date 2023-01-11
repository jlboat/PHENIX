\name{MPMM_helpers}
\alias{MPMM_helpers}
\alias{MPMM_Mstep}
\alias{MPMM_Estep}
\alias{MPMM_like}

\title{MPMM helper functions}
\description{The E step, M step and likelihood of the EM algorithm to fit an MPMM and plotting functions}

\usage{
MPMM_Mstep( Omega, pars, warm, init )
MPMM_Estep( C, D, Lam.R, Y.prime )
MPMM_like( Y.prime, Lam.R, C, D )
}

\arguments{
  \item{Omega}{An approximate sample covaiance matrix from E-step}
  \item{pars}{List of penalty parameters; type controls penalty type, rho controls penaltry extent}
  \item{warm}{Uses warm starts if true}
  \item{init}{Warm initialization}
  \item{C, D}{Genetic and environmental precision matrices}
  \item{Lam.R}{Eigenvalues of inverse kinship matrix}
  \item{Y.prime}{Kinship pre-whitened phenotype matrix}
}
\value{ 
'MPMM_Mstep' returns par, the solution to the penalized likelihood problem defined by the sufficient statistic Omega and input 'pars'. 'MPMM_Estep' returns C and D, estimates of the respective sample covariance, or 'Omega', matrices. 'MPMM_like' returns the likelihood of C and D given data Y.prime and Lam.R. 
}
\details{
'MPMM_Mstep' solves the optimization problem defined in equation (9) of the arxiv paper; this depends on the penalty function used, which is captured by 'pars'. 'MPMM_Estep' evaluates the expected sufficient statistics, given in expressions (7) and (8) in the arxiv paper. 'MPMM_like' evaluates the (unpenalized) likelihood of any 'C' and 'D' pair.
}

\author{Andy Dahl}

\keyword{ algebra }

\references{ Andy Dahl, Victoria Hore, Valentina Iotchkova, Jonathan Marchini (2013). Network inference in matrix-variate Gaussian models with non-independent noise. \url{http://arxiv.org/pdf/1312.1622.pdf} }

\seealso{\code{\link{MPMM}}, \code{\link{MPMM_impute}}}
