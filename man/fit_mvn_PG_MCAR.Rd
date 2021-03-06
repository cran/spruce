% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_mvn_PG_MCAR.R
\name{fit_mvn_PG_MCAR}
\alias{fit_mvn_PG_MCAR}
\title{Multivariate normal spatial mixture model clustering w/ PG multinomial regression on membership probabilities}
\usage{
fit_mvn_PG_MCAR(
  Y,
  W,
  coords_df,
  K,
  nsim = 2000,
  burn = 1000,
  z_init = NULL,
  verbose = FALSE
)
}
\arguments{
\item{Y}{An n x g matrix of gene expression values. n is the number of cell spots and g is the number of features.}

\item{W}{An n x v matrix of covariates to predict cluster membership. Should include an intercept (i.e., first column is 1)}

\item{coords_df}{An n x 2 data frame or matrix of 2d spot coordinates.}

\item{K}{The number of mixture components to fit.}

\item{nsim}{Number of total MCMC iterations to run.}

\item{burn}{Number of MCMC iterations to discard as burn in. The number of saved samples is nsim - burn.}

\item{z_init}{Optional initialized allocation vector. Randomly initialized if NULL.}

\item{verbose}{Logical for printing cluster allocations at each iteration.}
}
\value{
a list of posterior samples
}
\description{
Implement Gibbs sampling for MVN model with MCAR spatial random effects w/ PG multinomial regression on membership probabilities
}
\examples{
\donttest{
# parameters
data(coords_df_sim)
coords_df <- coords_df_sim[,1:2]
z <- remap_canonical2(coords_df_sim$z)
                                 
n <- nrow(coords_df) # number of observations
g <- 3 # number of features
K <- length(unique(coords_df_sim$z)) # number of clusters (mixture components)
pi <- table(z)/length(z) # cluster membership probability

W <- matrix(0, nrow = n, ncol = 2)
W[,1] <- 1
W[,2] <- sample(c(0,1),size = n, replace = TRUE, prob = c(0.5,0.5))

# Cluster Specific Parameters
Mu <- list(
  Mu1 = rnorm(g,-5,1),
  Mu2 = rnorm(g,0,1),
  Mu3 = rnorm(g,5,1),
  Mu4 = rnorm(g,-2,3)
)
# cluster specific variance-covariance
S <- matrix(1,nrow = g,ncol = g) # y covariance matrix
diag(S) <- 1.5
Sig <- list(
  Sig1 = S,
  Sig2 = S, 
  Sig3 = S,
  Sig4 = S
)

Y <- matrix(0, nrow = n, ncol = g)
for(i in 1:n)
{
  Y[i,] <- mvtnorm::rmvnorm(1,mean = Mu[[z[i]]],sigma = Sig[[z[i]]])
}

# fit model
# in practice use more mcmc iterations
fit <- fit_mvn_PG_MCAR(Y = Y, coords_df = coords_df, W = W, K = K, nsim = 10, burn = 0)}
}
