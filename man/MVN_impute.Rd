\name{MVN_impute}
\alias{MVN_impute}

\title{Matrix imputation with multivariate normal rows}
\description{An EM algorithm to iteratively impute missing values and estimate the column mean and covariance.}
\usage{
MVN_impute( Y, reltol=1e-4, intercept = TRUE, maxit=1e2, trace=FALSE )
}
\arguments{
  \item{Y}{A partially observed matrix}
  \item{reltol}{The algorithm terminates when the relative mean squared difference between successive imputations is less than reltol}
  \item{intercept}{Whether to fit column means as well as covariances}
  \item{maxit}{Maximum number of EM iterations}
  \item{trace}{Whether to print updates after each EM iteration}
}
\value{
  \item{Y}{The imputed matrix}
  \item{Sigma}{The estimated covariance across columns of Y}
}
\details{This models rows of Y as i.i.d. N( mu, Sigma ). The EM algorithm iteratively predicts the sufficient statistics and maximizes the resulting likelihood, ultimately returning imputed values for missing entries in Y and MLEs for mu and Sigma. If 'intercept'=FALSE, mu is constrained to 0.}

\examples{

  # simulate data
  N   <- 1e3
  P   <- 10
  Y0  <- matrix( rnorm( N*P ), N, P )
  Y   <- rmask( Y0, miss=.05 )

  # impute
  out <- MVN_impute( Y, trace=TRUE )
  head( out$Y )

}
\author{ Andy Dahl }

\references{ Andrew Dahl,  Valentina Iotchkova,  Amelie Baud,  Asa Johansson,  Ulf Gyllensten, Nicole Soranzo, Richard Mott, Andreas Kranis \& Jonathan Marchini (2016) A multiple-phenotype imputation method for genetic studies, \url{http://www.nature.com/ng/journal/v48/n4/full/ng.3513.html}, Nature Genetics, 48(4), 466-472 }
