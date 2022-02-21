#' Plot delta parameters from multinomial regression model
#'
#' Allows for visualization of multinomial regression models from spatial or non-spatial models
#'
#' @param fit A model fit returned by one of the fit_*_PG model functions
#'
#' @return a ggplot
#' @export
#' @import patchwork
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer separate
#' @importFrom dplyr filter
#' @importFrom tidyselect everything
#' @import ggplot2
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
#' plot_deltas(fit)}
plot_deltas <- function(fit)
{
  DELTA <- fit$DELTA
  v <- ncol(fit$W)
  K <- fit$K
  # name columns
  c_count = 1
  c_names <- rep("C",ncol(DELTA))
  for(k in 2:K)
  {
    for(vs in 1:v)
    {
      c_names[c_count] <- paste0("k",k,"_p",vs)
      c_count = c_count + 1
    }
  }
  colnames(DELTA) <- c_names
  DELTA <- as.data.frame(DELTA) 
  D_df_long <- tidyr::pivot_longer(DELTA,cols = tidyselect::everything(),names_to = "delta",values_to = "draw") 
  D_df_long <- tidyr::separate(data = D_df_long,col = "delta",into = c("k","p"),remove = FALSE)
  
  D_df_long_p1 <- dplyr::filter(D_df_long,.data$p == "p1")
  D1_add <- data.frame(delta = rep("k1_p1",nrow(DELTA)),
                       k = rep("k1",nrow(DELTA)),
                       p = rep("p1",nrow(DELTA)),
                       draw = rep(0,nrow(DELTA)))
  D_df_long_p1 <- rbind(D_df_long_p1,D1_add)
  d_plot_p <- ggplot(data = D_df_long_p1, aes(x = k, y = .data$draw, fill = .data$delta)) +
    geom_boxplot() + 
    geom_hline(yintercept = 0,linetype = "dashed") + 
    theme_classic() + 
    ylab("p = 1") + 
    xlab(NULL) + 
    annotate("text", x = 1,y = 0.1,label = "Reference") + 
    theme(legend.position = "none")
  d_plots <- d_plot_p
  for(vs in 2:v)
  {
    D_df_long_p <- dplyr::filter(D_df_long,.data$p == paste0("p",vs))
    D_add <- data.frame(delta = rep("k1_p1",nrow(DELTA)),
                         k = rep("k1",nrow(DELTA)),
                         p = rep(paste0("p",vs),nrow(DELTA)),
                         draw = rep(0,nrow(DELTA)))
    D_df_long_p <- rbind(D_df_long_p,D_add)
    d_plot_p <- ggplot(data = D_df_long_p, aes(x = k, y = .data$draw, fill = .data$delta)) +
      geom_boxplot() + 
      geom_hline(yintercept = 0,linetype = "dashed")+ 
      theme_classic() + 
      ylab(paste0("p = ",vs)) + 
      xlab(NULL)+ 
      theme(legend.position = "none")
    if(vs == v) d_plot_p <- d_plot_p + xlab("Cluster")
    d_plots <- d_plots + d_plot_p
  }
  d_plots <- d_plots + patchwork::plot_layout(nrow = v)
  return(d_plots)
}