\name{phenix_vb}

\alias{phenix_vb}
\alias{phenix_likelihood}

\title{Variational Bayes algorithm to fit phenix}
\description{Workhorse functions for phenix}

\usage{
phenix_vb( Y, Q, lam_K, N=nrow(Y), P=ncol(Y), M=min(c(N,P)),
    tau=0, e=P+5, E.inv=solve( diag(P)/(e-P-1) ),
    trace=0, reltol=1e-8, maxit=1e3, cutoff=1e-3, ... )
phenix_likelihood( Y, Q, lam_K, N, P, M, tau, E.inv,
    mu_Y, Sigma_Y, Omega, mu_beta, mu_S, eprime,
    svd.Omega, svd.V_S, svd.V_beta, svd.Omega_beta,
    j2miss, j2len, n.miss.types )
}

\arguments{
  \item{Y}{A matrix of phenotypes}
  \item{Q}{Eigenvectors of the kinship matrix}
  \item{lam_K}{Eigenvalues of the kinship matrix}
  \item{N}{Number of samples, or rows of Y}
  \item{P}{Number of phenotypes, or columns of Y}
  \item{M}{Maximum allowed rank of the low-rank, genetic term}
  \item{tau}{The prior precision on beta. Increasing tau increases shrinkage of the low-rank, genetic term}
  \item{e}{Weight of Wishart prior on environmental covariance; the prior is scaled to always have expectation diag('P')}
  \item{E.inv}{Inverse shape parameter of Wishart prior on environmental covariance}
  \item{trace}{Whether to print output at each VB iteration}
  \item{reltol}{Convergence criterion for relative change in the approximate log-likelihood}
  \item{maxit}{Maximum number of VB iterations}
  \item{cutoff}{Cutoff for component to be deemed zero}
  \item{mu_Y, Sigma_Y, Omega, mu_beta, mu_S, eprime}{Variational parameters}
  \item{svd.Omega, svd.V_S, svd.V_beta, svd.Omega_beta}{SVD's of variational parameters}
  \item{j2miss, j2len, n.miss.types}{Compressed representation of missingness pattern in Y}
  \item{...}{Arguments passed to 'MVN_impute', used for initialization}
}

\value{ 'phenix_vb' returns a locally optimal set of variational parameters for the approximation to the phenix posterior. 'phenix_likelihood' evaluates the marginal likelihood lower bound given by a set of variational parameters.}

\author{Andy Dahl}

\note{ These functions should rarely be called directly; the function 'phenix' should be used instead. }

\details{ 

'phenix_vb' is a deterministic algorithm that iteratively updates blocks of the variational parameters for the approximate posterior. The chosen variational approximation admits analytic solutions for each of these block updates; the algorithm iterates these updates (a.k.a. iterated conditional modes) until the relative change in the marginal likelihood lower bound (a.k.a. 'phenix_likelihood') is below 'reltol'.

'phenix_likelihood' is the objective function for 'phenix_vb'; maximizing this likelihood lower bound is equivalent to minimizing the KL-divergence between the exact posterior and its variational approximation. The VB algorithm is guaranteed to increase this likelihood at each iteration, and thus returns a local maximum of 'phenix_likelihood'.

}

\seealso{\code{\link{phenix}}}

\references{ Andrew Dahl,  Valentina Iotchkova,  Amelie Baud,  Asa Johansson,  Ulf Gyllensten, Nicole Soranzo, Richard Mott, Andreas Kranis \& Jonathan Marchini (2016) A multiple-phenotype imputation method for genetic studies, \url{http://www.nature.com/ng/journal/v48/n4/full/ng.3513.html}, Nature Genetics, 48(4), 466-472 }
