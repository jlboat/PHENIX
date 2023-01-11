\name{phenix}
\alias{phenix}
\title{PHENotype Imputation eXpedited}
\description{Imputes missing entries of a phenotype matrix and predicts breeding values}
\usage{
phenix( Y, K, Q, lam_K,
    test=FALSE, test.frac=.05, test.cols=1:ncol(Y),
    seed=8473, quantnorm=FALSE, scale=TRUE,
    trim=FALSE, trim.sds=4, ... )
}
\arguments{
  \item{Y}{A partially observed matrix of phenotypes}
  \item{K}{A positive semidefinite kinship matrix describing the genetic covariance between individuals}
  \item{Q}{Optional eigenvectors of the kinship matrix}
  \item{lam_K}{Optional eigenvalues of the kinship matrix}
  \item{test}{Option to test phenix's performance on a specific dataset by hiding some entries of Y and assessing recovery accuracy}
  \item{test.frac}{The fraction of Y's observed entries hidden in the testing procedure}
  \item{test.cols}{Which columns of Y can be masked in the testing procedure}
  \item{seed}{Enables reproducibility across test runs}
  \item{quantnorm}{Quantile normalizes phenotypes prior to imputation}
  \item{scale}{Whether to center and standardize columns of Y; this process is inverted after imputation so imputed phenotypes are on the input scale}
  \item{trim}{Whether to remove outliers, defined column-wise as entries more than 'trim.sds' standard deviations from their means}
  \item{trim.sds}{Number of standard deviations used to define outliers}
  \item{\dots}{Arguments passed to phenix_vb, mainly 'tau', which allows additional shrinkage of U, and 'maxit' and 'reltol', which define convergence}
}
\details{ phenix incorporates genetic relatedness into a Bayesian matrix factorization model to impute missing phenotypes in Y. It computes a low-rank reconstruction of U, the genetic contribution to phenotypic variance. It also estimates B and E, genetic and environmental phenotypic correlation matrices, from U and the environmental residual Y-U, respectively. }
\value{
  \item{imp}{Imputed phenotype matrix}
  \item{U}{Matrix of estimated breeding values}
  \item{S}{Matrix of latent components}
  \item{h2}{Vector of estimated phenotype heritabilities}
  \item{B}{Estimated genetic covariance matrix between phenotypes}
  \item{E}{Estimated environmental covariance matrix between phenotypes}
  \item{Q}{Eigenvectors of K}
  \item{lam_K}{Eigenvalues of K}
  \item{vb_out}{Output from phenix_vb, the workhorse function}
  \item{outliers}{The entries of 'c(Y)' set to NA prior to imputation}
  \item{test}{If test=TRUE, a list containing: 'time', the runtime; 'cors' and 'cor_glob', column-wise and overall imputation correlations; and 'mses' and 'mse_glob', column-wise and overall imputation mean-squared errors}
}
%\references{}
\author{Andy Dahl}
\examples{

  # simulate data
  N   <- 500
  P   <- 20
  K   <- sib_K( N=N )               # indepedendent sets of sibs
  Y0  <- rpheno( N=N, P=P, K=K )    # simulate phenotype matrix
  Y   <- rmask( Y0, miss=.05 )      # hide phenotypes

  out     <- phenix(Y,K,test=TRUE)  # impute missing data
  sub   <- which( is.na(Y) )
  cor( out$imp[ sub ], Y0[ sub ] )  #true overall accuracy
  out$test$cor_glob                 #estimated accuracies

}
\seealso{\code{\link{sim_data}},\code{\link{phenix_vb}}}

\references{ Andrew Dahl,  Valentina Iotchkova,  Amelie Baud,  Asa Johansson,  Ulf Gyllensten, Nicole Soranzo, Richard Mott, Andreas Kranis \& Jonathan Marchini (2016) A multiple-phenotype imputation method for genetic studies, \url{http://www.nature.com/ng/journal/v48/n4/full/ng.3513.html}, Nature Genetics, 48(4), 466-472 }
