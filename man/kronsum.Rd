\name{kronsum}

\alias{kronsum}
\alias{KS_dot_vec}
\alias{trp_KS}
\alias{tr}

\title{Kronecker Sum Operations with Low Complexity}
\description{Efficiently perform partial trace and matrix-vector multiplications with a (possibly inverse) Kronecker sum of symmetric matrices}

\usage{
KS_dot_vec( A, B, C, inv = FALSE, svd.A = svd(A), svd.C = svd(C), diag.A=FALSE )
trp_KS( keep, kill, svd.keep = svd(keep), lam.kill = svd(kill)$d ) 
tr(X)
}

\arguments{
  \item{A}{The N x N, symmetric left Kronecker summand }
  \item{C}{The P x P, symmetric right Kronecker summand }
  \item{B}{An N x P matrix }
  \item{X}{A square matrix }
  \item{inv}{Whether the Kronecker sum is inverted}
  \item{svd.A}{Optional svd of A }
  \item{svd.C}{Optional svd of C }
  \item{diag.A}{Whether A is diagonal}
  \item{keep}{The Kronecker summand whose dimension is returned by the partial trace }
  \item{kill}{The Kronecker summand whose dimension is traced out }
  \item{svd.keep}{The optional svd of the kept Kronecker summand }
  \item{lam.kill}{The eigenvalues of the killed Kronecker summand }
}

\value{ 'KS_dot_vec' returns a matrix of dimension dim(B); 'trp_KS' returns a matrix of dimension dim(keep); 'tr' returns a real number. }

\author{Andy Dahl}

\references{ Andy Dahl, Victoria Hore, Valentina Iotchkova, Jonathan Marchini (2013). Network inference in matrix-variate Gaussian models with non-independent noise. \url{http://arxiv.org/pdf/1312.1622.pdf} }
\note{ There are two partial traces associated with a Kronecker sum--this function computes the partial trace that has the same dimension as the 'keep' matrix.}

\details{

'KS_dot_vec' computes 

matrix( ( C \\oplus A )^i vec( B ) )

where i is -1 if inv=TRUE and 1 otherwise. The complexity is O( ( dim(A) + dim(C) )^3 ). If the either svd (of C or A) is provided, the respective cubic term(s) drop out. If, instead, A is diagonal, the complexity simply becomes O( dim(A) dim(C)^2 + dim(C)^3 ).

'trp_KS' computes

trp( ( keep \\oplus kill )^i )

where i is as above. Here, trp is the partial trace that returns a matrix of the same dimensions as 'keep'. The complexity is O( dim(keep) dim(kill) + dim( keep )^3 ) if the eigendecomposition of 'keep' is given.

'tr' computes the trace of a matrix.

}

\examples{
  A <- rWishart( 1, 15, diag(7) )[,,1]
  B <- matrix( rnorm( 7*11 ), 7, 11 )
  C <- rWishart( 1, 15, diag(11) )[,,1]

  C_oplus_A <- ( C \%x\% diag(7) + diag(11) \%x\% A )

  # compare KS_dot_vec to naive implementation
  res1  <- KS_dot_vec( A, B, C, inv=TRUE )
  res2  <- matrix( solve( C_oplus_A ) \%*\% c( B ), 7, 11 )
  print( all.equal( res1, res2 ) )

  # compare trp_KS to naive implementation (on entry [3,4] )
  res1  <- trp_KS( C, A )
  rows  <- 1:nrow(A) + (3-1)*nrow(A)
  cols  <- 1:nrow(A) + (4-1)*nrow(A)
  res2  <- tr( solve( C_oplus_A )[ rows, cols ] )
  print( all.equal( res1[3,4], res2 ) )
}

\keyword{algebra}
