#' Calculate cluster uncertainty
#'
#' Use posterior estimates to calculate uncertainty scores
#'
#' @param fit A model fit returned by one of the fit_*_PG model functions
#'
#' @return An n x (K + 1) matrix. First K columns are continuous phenotypes, and last column is uncertainty scores 
#' @export
#' @importFrom mvtnorm dmvnorm
#' @examples 
#' \donttest{
#' # parameters
#' data(coords_df_sim)
#' coords_df <- coords_df_sim[,1:2]
#' z <- remap_canonical2(coords_df_sim$z)
#'                                  
#' n <- nrow(coords_df) # number of observations
#' g <- 3 # number of features
#' K <- length(unique(coords_df_sim$z)) # number of clusters (mixture components)
#' pi <- table(z)/length(z) # cluster membership probability
#' 
#' W <- matrix(0, nrow = n, ncol = 2)
#' W[,1] <- 1
#' W[,2] <- sample(c(0,1),size = n, replace = TRUE, prob = c(0.5,0.5))
#' 
#' # Cluster Specific Parameters
# cluster specific means
#' Mu <- list(
#'   Mu1 = rnorm(g,-5,1),
#'   Mu2 = rnorm(g,0,1),
#'   Mu3 = rnorm(g,5,1),
#'   Mu4 = rnorm(g,-2,3)
#' )
#' # cluster specific variance-covariance
#' S <- matrix(1,nrow = g,ncol = g) # y covariance matrix
#' diag(S) <- 1.5
#' Sig <- list(
#'   Sig1 = S,
#'   Sig2 = S, 
#'   Sig3 = S,
#'   Sig4 = S
#' )
#' 
#' Y <- matrix(0, nrow = n, ncol = g)
#' for(i in 1:n)
#' {
#'   Y[i,] <- mvtnorm::rmvnorm(1,mean = Mu[[z[i]]],sigma = Sig[[z[i]]])
#' }
#' 
#' # fit model
#' # in practice use more mcmc iterations
#' fit <- fit_mvn_PG_smooth(Y = Y, coords_df = coords_df, W = W, K = K, nsim = 10, burn = 0)
#' scores_df <- get_scores(fit)}
#' 
#' 
get_scores <- function(fit)
{
  Y <- fit$Y
  n <- nrow(Y)
  W <- fit$W
  v <- ncol(W)
  MU <- fit$MU
  K <- fit$K
  p <- ncol(Y)
  SIGMA <- fit$SIGMA
  DELTA <- fit$DELTA
  z_map <- fit$z
  
  delta_post <- matrix(colMeans(DELTA),nrow = v,ncol = K-1)
  eta_post <- cbind(rep(0,n),W%*%delta_post)
  PI_post <- exp(eta_post)/(1+apply(as.matrix(exp(eta_post[,-1])),1,sum))
  pi_post <- table(z_map)/n
  PI_d_post <- matrix(0,nrow = n,ncol = K)
  
  mu_post <- lapply(MU, colMeans)
  sig_post <- lapply(SIGMA, colMeans)
  
  ro <- matrix(0,nrow = n,ncol = K)
  
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  for(i in 1:n)
  {
    pj <- rep(0,K)
    yi <- Y[i,]
    for(k in 1:K)
    {
      mu_k_post <- mu_post[[k]]
      sigma_k_post <- matrix(sig_post[[k]],
                             nrow = p,
                             ncol = p)
      pj[k] <- dmvnorm(yi,mu_k_post,sigma_k_post)
    }
    pii <- PI_post[i,]
    pjn <- pii*pj/sum(pii*pj)
    PI_d_post[i,] <- pjn
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  u_score <- rep(0,n)
  for(i in 1:n)
  {
    u_score[i] <- 1-max(PI_d_post[i,])
  }
  
  ret <- cbind(PI_d_post,u_score)
  colnames(ret) <- c(paste0("k",1:K,"_score"),"u_score")
  
  return(ret)
}