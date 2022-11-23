#'
#' This function generates a sample from the posterior of COMIX.
#'
#' @param Y Matrix of the data. Each row represents an observation.
#' @param C Vector of the group label of each observation. Labels must be integers starting from 1.
#' @param prior A list giving the prior information. If unspecified, a default prior is used.
#' The list includes the following parameters:
#' \itemize{
#' \item \code{zeta}: Coarsening parameter. A number between 0 and 1. \code{zeta} = 1: sample from standard posterior;
#' \code{zeta} < 1: sample from power posterior. The lower \code{zeta} is, the more flexible the kernels become.
#' \item \code{K}: Maximal number of mixture components.
#' \item \code{eta_prior} Parameters for gamma prior for concentration parameter of the stick breaking process
#' prior for the weights.
#' \item \code{m0}: Number of degrees of freedom for the inverse Wishart prior for Sigma, the covariance matrix
#' of the kernels. 
#' \item \code{Lambda}: Mean parameter for the inverse Wishart prior for Sigma, the covariance matrix
#' of the kernels. 
#' \item \code{b0}: Mean parameter for the multivariate normal distribution that is the prior for the
#' group mean parameter xi0.
#' \item \code{B0}: Covariance parameter for the multivariate normal distribution that is the prior for the
#' group mean parameter xi0.
#' \item \code{e0}: Number of degrees of freedom for the inverse Wishart prior for \eqn{E_k}{E_k}, the 
#' covariance matrix of the multivariate normal from which \eqn{\xi_{j,k}}{xi_jk} are drawn.
#' \item \code{E0}: Mean parameter for the inverse Wishart prior for \eqn{E_k}{E_k}, the 
#' covariance matrix of the multivariate normal from which \eqn{\xi_{j,k}}{xi_jk} are drawn.
#' \item \code{merge_step}: Introduce step to merge mixture components with small KL divergence. Default 
#' is \code{merge_step = TRUE}.
#' \item \code{merge_par}: Parameter controlling merging radius. Default is \code{merge_par = 0.1}.
#' }
#' @param pmc A list giving the Population Monte Carlo (PMC) parameters:
#' \itemize{
#' \item \code{npart}: Number of PMC particles.
#' \item \code{nburn}: Number of burn-in steps
#' \item \code{nsave}: Number of steps in the chain after burn-in.
#' \item \code{nskip}: Thinning parameter, number of steps to skip between saving steps after burn-in.
#' \item \code{ndisplay}: Display status of chain after every \code{ndisplay} steps.
#' }
#' @param state A list giving the initial cluster labels:
#' \itemize{
#' \item \code{t}: An integer vector, same length as the number of rows of \code{Y}, with cluster labels
#' between \code{1} and \code{K}.
#' }
#' @param ncores The number of CPU cores to utilize in parallel. Defaults to 2.
#' @return An object of class COMIX, a list of 4:
#' \itemize{
#' \code{chain}, a named list:
#' \itemize{
#' \item  \code{t}: an \code{nsave} \eqn{\times}{x} \code{nrow(Y)} matrix with estimated cluster labels
#' for each saved step of the chain and each observation in the data \code{Y}.
#' \item  \code{z}: a \code{nsave} \eqn{\times}{x} \code{nrow(Y)} matrix with estimated values of 
#' the latent \eqn{z_{i,j}}{z_ij} variable for the parameterization of the
#' multivariate skew normal distribution used in the sampler for each saved step of 
#' the chain and each observation in the data \code{Y}.
#' \item  \code{W}: an \code{length(unique(C))} \eqn{\times}{x} \code{K} \eqn{\times}{x} 
#' \item  \code{nsave}: array storing the estimated sample- and cluster-specific weights for each 
#' saved step of the chain.
#' \item  \code{xi}: an \code{length(unique(C))} \eqn{\times}{x} \code{(ncol(Y) x K)} 
#' \eqn{\times}{x} \code{nsave} array storing the estimated sample- and cluster-specific
#' multivariate skew normal location parameters of the kernel for each saved step of the chain.
#' \item  \code{xi0}: an \code{ncol(Y)} \eqn{\times}{x} \code{K} \eqn{\times}{x} 
#' \item  \code{nsave}: array storing the estimated cluster-specific 
#' group location parameters for each saved step of the chain.
#' \item  \code{psi}: an \code{ncol(Y)} \eqn{\times}{x} \code{K} \eqn{\times}{x} \code{nsave}
#' array storing the estimated cluster-specific skew parameters of the kernels in
#' the parameterization of the
#' multivariate skew normal distribution used in the sampler
#' for each saved step of the chain.
#' \item  \code{G}: an  \code{ncol(Y)} \eqn{\times}{x} \code{(ncol(Y) x K)} 
#' \eqn{\times}{x} \code{nsave} array storing the estimated cluster-specific
#' multivariate skew normal scale matrix (in row format) of the kernel
#' used in the sampler for each saved step of the chain.
#' \item  \code{E}: an  \code{ncol(Y)} \eqn{\times}{x} \code{(ncol(Y) x K)} 
#' \eqn{\times}{x} \code{nsave} array storing the estimated covariance matrix 
#' (in row format) of the multivariate normal distribution from which the  
#' sample- and cluster-specific location parameters are drawn for each saved step 
#' of the chain.
#' \item \code{eta}: a \code{nsave} \eqn{\times}{x} \code{1} matrix storing the
#' estimated Dirichlet Process Mixture concentration parameter for each saved step of the chain.
#' \item \code{Sigma}: an  \code{ncol(Y)} \eqn{\times}{x} \code{(ncol(Y) x K)} 
#' \eqn{\times}{x} \code{nsave} array storing the estimated cluster-specific
#' multivariate skew normal scale matrix (in row format) of the kernel for each saved step of the chain.
#' \item \code{alpha}: an \code{ncol(Y)} \eqn{\times}{x} \code{K} \eqn{\times}{x} \code{nsave}
#' array storing the estimated cluster-specific skew parameters of the
#' kernel's multivariate skew normal distribution
#' for each saved step of the chain.
#' }
#' \item \code{data}, a named list that includes the matrix of the data \code{Y}
#' and \code{C} the vector of the group label of each observation.
#' \item \code{prior} and \code{pmc}, the lists, as above, that were provided as inputs to 
#' the function.
#' }
#' @examples
#' library(COMIX)
#' # Number of observations for each sample (row) and cluster (column):
#' njk <- 
#'   matrix(
#'     c(
#'       150, 300,
#'       250, 200
#'     ),
#'     nrow = 2,
#'     byrow = TRUE
#'   )
#' 
#' 
#' # Dimension of data:
#' p <- 3
#' 
#' # Scale and skew parameters for first cluster:
#' Sigma1 <- matrix(0.5, nrow = p, ncol = p) + diag(0.5, nrow = p)
#' alpha1 <- rep(0, p)
#' alpha1[1] <- -5
#' # location parameter for first cluster in first sample:
#' xi11 <- rep(0, p)
#' # location parameter for first cluster in second sample (aligned with first):
#' xi21 <- rep(0, p)
#' 
#' # Scale and skew parameters for second cluster:
#' Sigma2 <- matrix(-1/3, nrow = p, ncol = p) + diag(1 + 1/3, nrow = p)
#' alpha2 <- rep(0, p)
#' alpha2[2] <- 5
#' # location parameter for second cluster in first sample:
#' xi12 <- rep(3, p)
#' # location parameter for second cluster in second sample (misaligned with first):
#' xi22 <- rep(4, p)
#' 
#' # Sample data:
#' set.seed(1)
#' Y <- 
#'   rbind(
#'     sn::rmsn(njk[1, 1], xi = xi11, Omega = Sigma1, alpha = alpha1),
#'     sn::rmsn(njk[1, 2], xi = xi12, Omega = Sigma2, alpha = alpha2),
#'     sn::rmsn(njk[2, 1], xi = xi21, Omega = Sigma1, alpha = alpha1),
#'     sn::rmsn(njk[2, 2], xi = xi22, Omega = Sigma2, alpha = alpha2)
#'   )
#' 
#' C <- c(rep(1, rowSums(njk)[1]), rep(2, rowSums(njk)[2]))
#' 
#' prior <- list(zeta = 1, K = 10)
#' pmc <- list(naprt = 5, nburn = 200, nsave = 200) # Reasonable usage
#' pmc <- list(naprt = 5, nburn = 2, nsave = 5) # Minimal usage for documentation
#' # Fit the model:
#' res <- comix(Y, C, pmc = pmc, prior = prior)
#' 
#' # Relabel to resolve potential label switching issues:
#' res_relab <- relabelChain(res)
#' 
#' # Generate calibrated data:
#' cal <- calibrateNoDist(res_relab)
#' 
#' # Compare raw and calibrated data: (see plot in vignette)
#' # par(mfrow=c(1, 2))
#' # plot(Y, col = C, xlim = range(Y[,1]), ylim = range(Y[,2]) )
#' 
#' # Get posterior estimates for the model parameters:
#' res_summary <- summarizeChain(res_relab)
#' # Check for instance, the cluster assignment labels:
#' table(res_summary$t)
#' # Indeed the same as 
#' colSums(njk)
#' 
#' # Or examine the skewness parameter for the non-trivial clusters:
#' res_summary$alpha[ , unique(res_summary$t)]
#' # And compare those to
#' cbind(alpha1, alpha2)
#' 
#' # (see vignette for a more detailed example)
#' @export
comix = function(Y, C, prior = NULL, pmc = NULL, state = NULL, ncores = 2)
{
  Y = as.matrix(Y)
  p = ncol(Y)
  
  # R wrapper:
  if(is.null(prior)) {
    prior = list(    zeta = 1,
                     K = 10,
                     eta_prior = c(1, 1),
                     m0 = ncol(Y) + 2,
                     Lambda = stats::cov(Y),
                     b0 = colMeans(Y),
                     B0 = 100 * stats::cov(Y),
                     e0 = ncol(Y) + 2,
                     E0 = 0.1 * stats::cov(Y),
                     merge_step = TRUE,
                     merge_par = 0.1)
  } else {
    if(is.null(prior$zeta)) 
      prior$zeta = 1;
    if(is.null(prior$K)) 
      prior$K = 10;
    if(is.null(prior$eta_prior))
      prior$eta_prior = c(1, 1);
    if(is.null(prior$Lambda))
      prior$Lambda = stats::cov(Y);
    if(is.null(prior$m0))
      prior$m0 = ncol(Y) + 2;
    if(is.null(prior$b0))
      prior$b0 = colMeans(Y);
    if(is.null(prior$B0))
      prior$B0 = 100*stats::cov(Y);
    if(is.null(prior$e0))
      prior$e0 = ncol(Y) + 2;
    if(is.null(prior$E0))
      prior$E0 = 0.1 * stats::cov(Y);
    if(is.null(prior$merge_step))
      prior$merge_step = TRUE;
    if(is.null(prior$merge_par))
      prior$merge_par = 0.1;
  }
  
  
  if(is.null(pmc)) {
    pmc = list(npart = 10, nburn = 1000, nsave = 1000, nskip = 1, ndisplay = 500)
  } else {
    if(is.null(pmc$npart))
      pmc$npart = 10
    if(is.null(pmc$nburn))
      pmc$nburn = 5000
    if(is.null(pmc$nburn))
      pmc$nburn = 5000
    if(is.null(pmc$nsave))
      pmc$nsave = 1000
    if(is.null(pmc$nskip))
      pmc$nskip = 1
    if(is.null(pmc$ndisplay))
      pmc$ndisplay = 100
  }

  if(is.null(state$t)) {
    state$t = stats::kmeans(Y, prior$K, iter.max = 100)$cluster
  }

  J = length(unique(C))
  if( sum( sort(unique(C)) == 1:J )  != J )
  {
    print("ERROR: unique(C) should look like 1, 2, ...")
    return(0);
  }
  C = C - 1
  state$t = state$t - 1
  
  ans = perturbedSNcpp(Y, C, prior, pmc, state, initParticles = NULL, init=T, ncores=2)
  colnames(ans$data$Y) = colnames(Y)
  ans$data$C = ans$data$C + 1
  ans$chain$t = ans$chain$t + 1
  class(ans) = "COMIX"
  return(ans)
}
