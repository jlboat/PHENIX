\name{misc}

\alias{misc}
\alias{mat.sqrt}
\alias{schur}
\alias{vec.pinv}
\alias{quantnorm}
\alias{logdet}

\title{Miscellaneous support functions for phenix}
\description{Miscellaneous support functions for phenix}

\usage{
mat.sqrt( X, inv = FALSE )
schur( X, m )
vec.pinv( x, tol = 1e-08 )
quantnorm( Y )
logdet( X )
}

\arguments{
  \item{X}{A symmetric, positive semidefinite matrix}
  \item{inv}{Whether to return the inverted square root}
  \item{m}{An integer vector of indices defining the Schur complement submatrix}
  \item{x}{A numeric vector}
  \item{tol}{Value at which entries of x are deemed numerically 0}
  \item{Y}{A partially observed matrix}
}

\value{ 'mat.sqrt' returns a matrix of dimension dim(X). 'schur' returns the Schur complement of X[m,m] in X, an |m| x |m| matrix. 'vec.pinv' performs element-wise pseudoinversion. 'quantnorm' returns a matrix, of the same dimension and missingness pattern as Y, with quantile-normalized columns. 'logdet' returns a scalar. }

\author{Andy Dahl}

\details{ 

'mat.sqrt' square-roots a covariance matrix by square-rooting (and potentially inverting) its eigenvalues; if any eigenvalues are non-positive (X is not psd), they are set to 1e-8. 

'vec.pinv' performs element-wise pseudoinversion; this enables psd matrix pseudoinversion by applying 'vec.pinv' to eigenvalues. 

'quantnorm' matches the quantiles of (the observed entries of) columns of Y to standard normal and then scales them to be mean zero and unit variance (which does not automatically hold in the presence of ties). 

}

\examples{

  N     <- 6
  P     <- 2
  (K    <- sib_K(N,fam_size=3))

  mat.sqrt(K)
  schur(K, 1:5)
  logdet(K)

  (Y     <- rpheno(N,P,K))
  (Ymiss <- rmask(Y, miss= 0.2))
  vec.pinv(Y[,1], tol = 1e-08)
  quantnorm(Y)

}
