\name{MPMM}

\alias{MPMM}
\alias{lik.fxn}

\title{Fit a MultiPhenotype Mixed Model}

\usage{
MPMM( Y, K, Y.prime, Lam.R, tol=1e-4, maxit=1e4, warm.updates=TRUE,
       trace=FALSE, print.like=trace, 
       lmm_init=TRUE, start.C, start.D, eig.tol=1e-6, standardize=missing(Y.prime) )
}

\arguments{
  \item{Y}{Raw, fully observed phenotype matrix}
  \item{K}{A positive semidefinite kinship matrix describing the genetic covariance between individuals }
  \item{Y.prime}{Kinship pre-whitened phenotype matrix}
  \item{Lam.R}{Eigenvalues of inverse kinship matrix}
  \item{tol}{Convergence criterion}
  \item{maxit}{Maximum number of EM iterations}
  \item{warm.updates}{Initializes M steps with solutions from last iteration if true}
  \item{trace}{Prints output if true}
  \item{print.like}{Includes likelihood evaluations in printout}
  \item{lmm_init}{Initializes C and D with univariate mixed model}
  \item{start.C, start.D}{Optional initalizations}
  \item{eig.tol}{Eigenvalues of K below eig.tol are thresholded to 0}
  \item{standardize}{If TRUE and Y given and unscaled, a warning is produced and Y is scaled internally}
}
\value{
  \item{C}{(Locally) maximum likelihood estimate for genetic precision between phenotypes}
  \item{D}{(Locally) maximum likelihood estimate for environmental precision between phenotypes}
  \item{h2}{Estimated heritability from C and D}
  \item{times}{Cumulative times spent in E-step, M-step for C and M-step for D}
  \item{del}{Sequence of relative Frobenius changes in parameters}
  \item{niter}{Number of EM iterations}
  \item{conv}{Whether EM converged prior to 'maxit' iterations}
  \item{C.unscaled,D.unscaled}{C and D on the input Y scale of (if standardize=TRUE and Y is unscaled)} 
  \item{ll}{Log-likelihood of returned C and D}
  \item{ll_path}{Log-likelihood path through EM iterations (if requeested)}
}

\author{Andy Dahl}

\description{Fit local maximum likelihood variance components in a multiple phenotype mixed model.}

\details{ This function fits C and D using the model

Y ~ MN( 0, K, C^{-1} ) + MN( 0, I, D^{-1} )

where Y is a fully observed matrix of phenotypes, K is a known kinship matrix, and C and D describe the heritable and environmental correlations between phenotypes, respectively. There is also the option to add a variety of penalties to C and D: NA fits the unpenalized model; 'ind' requires the precision matrix to be diagonal; 'iid' requires it to be spherical; 'l1' and 'l2' add lasso ridge penalties, respectively, to all entries of the precision matrix; and 'l1_cor' adds a lasso penalty on the partial correlation matrix.

Y.prime is a version of Y that has been whitened by the eigenvectors of K (i.e. if K = U diag(1/Lam.R) U^T is an eigendecomposition, Y.prime = U^T Y), so that rows of Y.prime are independent.

The algorithm terminates after the relative Frobenius change in C and D is less than reltol or maxit iterations.

Essentially the only effect of the eigenvalue thresholding is to make liklihood computations more stable.

}

\note{The convexity of MPMM is not currently understood, and the output C and D are guaranteed only to be local maximizers of the likelihood.}

\examples{

  # simulate data
  N   <- 1e3
  P   <- 7
  h2  <- rep( .5, P ) # vector of heritabilities
  K   <- sib_K( N=N ) # kinship matrix of indepedendent sets of sibs
  Y   <- rpheno( N=N, P=P, K=K, h2=h2 ) # simulate phenotype matrix

  # fit an MPMM
  out     <- MPMM(scale(Y),K,trace=TRUE,maxit=1e3)
  out$h2

}
\keyword{ algebra }

\references{ Andy Dahl, Victoria Hore, Valentina Iotchkova, Jonathan Marchini (2013). Network inference in matrix-variate Gaussian models with non-independent noise. \url{http://arxiv.org/pdf/1312.1622.pdf} }
\seealso{\code{\link{MPMM_impute}}, \code{\link{MPMM_helpers}}}
