#' This function aligns multiple samples so that their location parameters
#' are equal.
#'
#' @param x An object of class COMIX.
#' @param reference.group An integer between 1 and the number of groups in the data
#' (\code{length(unique(C))}). Defaults to \code{NULL}. If \code{NULL}, the samples 
#' are aligned so that their location parameters are set to be at the estimated
#' group location parameter. If an integer, the samples are aligned so that their 
#' location parameters are the same as the location parameter of sample \code{reference.group}.
#' @return A named list of 3:
#' \itemize{
#' \item \code{Y_cal}: a \code{nrow(x$data$Y)} \eqn{\times}{x} \code{ncol(x$data$Y)}
#' matrix, a calibrated version of the original data.
#' \item \code{calibration_distribution}: an \code{x$pmc$nsave} \eqn{\times}{x} 
#' \code{ncol(x$data$Y)} \eqn{\times}{x} \code{nrow(x$data$Y)} array storing the 
#' difference between the estimated sample-specific location parameter and the group
#' location parameter for each saved step of the chain.
#' \item \code{calibration_median}: a \code{nrow(x$data$Y)} \eqn{\times}{x} \code{ncol(x$data$Y)}
#' matrix storing the median difference between the estimated sample-specific location parameter and the group
#' location parameter for each saved step of the chain. This matrix is equal to the 
#' difference between the uncalibrated data (\code{x$data$Y}) and the calibrated
#' data (\code{Y_cal}).
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
#' @export
calibrate <- function(x, reference.group = NULL)
{
  C = x$data$C - 1
  Z = x$chain$t - 1
  ref = ifelse(is.null(reference.group), -1, reference.group - 1 )
  
  if (!is.matrix(Z)) Z=matrix(Z, ncol=1)
  ns  = dim(x$chain$xi0)[3]
  K = dim(x$chain$xi0)[2]

  output = calib(x$data$Y,
                 matrix(C,ncol=1),
                 Z,
                 x$chain$xi, dim(x$chain$xi),
                 x$chain$xi0, dim(x$chain$xi0),
                 ref)
  colnames(output$Y_cal) = colnames(x$data$Y)
  return(output)
}

#' This function aligns multiple samples so that their location parameters
#' are equal.
#'
#' @param x An object of class COMIX.
#' @param reference.group An integer between 1 and the number of groups in the data
#' (\code{length(unique(C))}). Defaults to \code{NULL}. If \code{NULL}, the samples 
#' are aligned so that their location parameters are set to be at the estimated
#' group location parameter. If an integer, the samples are aligned so that their 
#' location parameters are the same as the location parameter of sample \code{reference.group}.
#' @return A named list of 2:
#' \itemize{
#' \item \code{Y_cal}: a \code{nrow(x$data$Y)} \eqn{\times}{x} \code{ncol(x$data$Y)}
#' matrix, a calibrated version of the original data.
#' \item \code{calibration_median}: a \code{nrow(x$data$Y)} \eqn{\times}{x} \code{ncol(x$data$Y)}
#' matrix storing the median difference between the estimated sample-specific location parameter and the group
#' location parameter for each saved step of the chain. This matrix is equal to the 
#' difference between the uncalibrated data (\code{x$data$Y}) and the calibrated
#' data (\code{Y_cal}).
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
#' @export
calibrateNoDist <- function(x, reference.group = NULL)
{
  C = x$data$C - 1
  Z = x$chain$t - 1
  ref = ifelse(is.null(reference.group), -1, reference.group - 1 )
  
  if (!is.matrix(Z)) Z=matrix(Z, ncol=1)
  ns  = dim(x$chain$xi0)[3]
  K = dim(x$chain$xi0)[2]
  
  output = calibNoDist(x$data$Y,
                 matrix(C,ncol=1),
                 Z,
                 x$chain$xi, dim(x$chain$xi),
                 x$chain$xi0, dim(x$chain$xi0),
                 ref)
  colnames(output$Y_cal) = colnames(x$data$Y)
  return(output)
  
}

#' This function relabels the chain to avoid label switching issues.
#'
#' @param res An object of class COMIX.
#' @return An object of class COMIX where \code{res$chain$t} is replaced with the
#' new labels.
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
#' @export
relabelChain = function(res) {
  res$chain$t = res$chain$t - 1
  relabeled_chain = relabel(res)
  res$chain = relabeled_chain 
  res
}

#' Convert between parameterizations of the multivariate skew normal distribution.
#'
#' @param Sigma A scale matrix.
#' @param alpha A vector for the skew parameter.
#' @return A list:
#' \itemize{
#' \item \code{delta}: a reparameterized skewness vector, a transformed 
#' version of \code{alpha}.
#' \item \code{omega}: a diagonal matrix of the same dimensions as \code{Sigma}, 
#' the diagonal elements are the square roots of the diagonal elements of \code{Sigma}.
#' \item \code{psi}: another reparameterized skewness vector, utilized in the sampler.
#' \item \code{G}: a reparameterized version of \code{Sigma}, utilized in the sampler.
#' }
#' @examples
#' library(COMIX)
#' # Scale and skew parameters:
#' Sigma <- matrix(0.5, nrow = 4, ncol = 4) + diag(0.5, nrow = 4)
#' alpha <- c(0, 0, 0, 5)
#' transformed_parameters <- transform_params(Sigma, alpha)
#' @export
transform_params = function(Sigma, alpha) {
  n = NROW(Sigma)
  m = NCOL(Sigma)
  if (n!=m) stop("Sigma is not sqaure")
  if (length(alpha)!=n) stop("alpha is of wrong length")
  if (n==1) {
    omega = sqrt(Sigma)
  } else {
    omega = sqrt(diag(diag(Sigma)))
  }
  omega_inv = solve(omega)
  Omega = omega_inv %*% Sigma %*% omega_inv
  alpha = matrix(alpha, ncol=1)
  numer = Omega %*% alpha
  denom = as.numeric(sqrt(1 + t(alpha) %*% Omega %*% alpha))
  delta = numer/denom
  if (n==1) {
    psi = sqrt(Sigma) * delta
  } else {
    psi = sqrt(diag(Sigma)) * delta
  }
  return(list(delta=delta, omega=omega,
              psi=psi, G = (Sigma - psi%*%t(psi))))
}


Mode = function(x) {
  ux = unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' This function provides post-hoc estimates of the model parameters.
#'
#' @param res An object of class COMIX.
#' @return A named list:
#' \itemize{
#' \item \code{xi0}: a \code{ncol(res$data$Y)} \eqn{\times}{x} \code{res$prior$K} matrix storing
#' the posterior mean of the group location parameter.
#' \item \code{psi}: a \code{ncol(res$data$Y)} \eqn{\times}{x} \code{res$prior$K} matrix storing
#' the posterior mean of the multivariate skew normal kernels skewness parameter (in the parameterization used in
#' the sampler).
#' \item \code{alpha}: a \code{ncol(res$data$Y)} \eqn{\times}{x} \code{res$prior$K} matrix storing
#' the posterior mean of the multivariate skew normal kernels skewness parameter.
#' \item \code{W}: a \code{length(unique(res$data$C))} \eqn{\times}{x} \code{res$prior$K} matrix storing
#' the posterior mean of the mixture weights for each sample and cluster.
#' \item \code{xi}: an \code{length(unique(res$data$C))} \eqn{\times}{x} \code{ncol(res$data$Y)}
#' \eqn{\times}{x} \code{res$prior$K} array storing the the posterior mean of the 
#' multivariate skew normal kernels location parameter for each sample and cluster.
#' \item \code{Sigma}: an \code{ncol(res$data$Y)} \eqn{\times}{x} \code{ncol(res$data$Y)}
#' \eqn{\times}{x} \code{res$prior$K} array storing the the posterior mean of the 
#' scaling matrix of the multivariate skew normal kernels for each cluster.
#' \item \code{G}: an \code{ncol(res$data$Y)} \eqn{\times}{x} \code{ncol(res$data$Y)}
#' \eqn{\times}{x} \code{res$prior$K} array storing the the posterior mean of the 
#' scaling matrix of the multivariate skew normal kernels for each cluster (in the 
#' parameterization used in the sampler).
#' \item \code{E}: an \code{ncol(res$data$Y)} \eqn{\times}{x} \code{ncol(res$data$Y)}
#' \eqn{\times}{x} \code{res$prior$K} array storing the the posterior mean of the 
#' covariance matrix of the multivariate normal distributions for each cluster form which
#' the sample specific location parameters are drawn.
#' \item \code{meanvec}: an \code{length(unique(res$data$C))} \eqn{\times}{x} \code{ncol(res$data$Y)}
#' \eqn{\times}{x} \code{res$prior$K} array storing the the posterior mean of the 
#' multivariate skew normal kernels mean parameter for each sample and cluster.
#' \item \code{meanvec0}: a \code{ncol(res$data$Y)} \eqn{\times}{x} \code{res$prior$K} matrix storing
#' the posterior mean of the group mean parameter.
#' \item \code{t}: Vector of length \code{nrow(x$data$Y)}. Each element is the mode
#' of the posterior distribution of cluster labels.
#' \item \code{eta}: scalar, the mean of the posterior distribution of the estimated
#' Dirichlet Process Mixture concentration parameter.
#'}
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
#' @export
summarizeChain = function( res ) {
  chainSummary = list()
  K = res$prior$K
  chain = res$chain
  p = ncol(res$data$Y)
  J = length(unique(res$data$C))
  ns = res$pmc$nsave
  
  xi_raw = matrix(0, nrow=J, ncol=p*K)
  chainSummary$xi0 = matrix(0, nrow=p,ncol=K)
  chainSummary$psi = matrix(0, nrow=p,ncol=K)
  Sigma_raw = matrix(0, nrow=p, ncol=p*K)
  G_raw = matrix(0, nrow=p, ncol=p*K)
  E_raw = matrix(0, nrow=p, ncol=p*K)
  chainSummary$alpha = matrix(0, nrow=p,ncol=K)
  chainSummary$W = matrix(0, nrow=J, ncol=K)
  for (i in 1:ns) {
    xi_raw = xi_raw + chain$xi[,,i]
    chainSummary$xi0 = chainSummary$xi0 + chain$xi0[,,i]
    chainSummary$psi = chainSummary$psi + chain$psi[,,i]
    Sigma_raw = Sigma_raw + chain$Sigma[,,i]
    G_raw = G_raw + chain$G[,,i]
    E_raw = E_raw + chain$E[,,i]
    chainSummary$alpha = chainSummary$alpha + chain$alpha[,,i]
    chainSummary$W = chainSummary$W + chain$W[,,i]
  }
  
  chainSummary$W = chainSummary$W/ns
  
  xi_raw = xi_raw/ns
  chainSummary$xi = array(0, dim=c(J,p,K))
  for (k in 1:K) {
    chainSummary$xi[,,k] = xi_raw[,(1+p*(k-1)):(p*k)]
  }
  
  chainSummary$xi0 = chainSummary$xi0/ns
  chainSummary$psi = chainSummary$psi/ns
  
  Sigma_raw = Sigma_raw/ns
  chainSummary$Sigma = array(0, dim=c(p,p,K))
  for (k in 1:K) {
    chainSummary$Sigma[,,k] = Sigma_raw[,(1+p*(k-1)):(p*k)]
  }
  
  G_raw = G_raw/ns
  chainSummary$G = array(0, dim=c(p,p,K))
  for (k in 1:K) {
    chainSummary$G[,,k] = G_raw[,(1+p*(k-1)):(p*k)]
  }
  
  E_raw = E_raw/ns
  chainSummary$E = array(0, dim=c(p,p,K))
  for (k in 1:K) {
    chainSummary$E[,,k] = E_raw[,(1+p*(k-1)):(p*k)]
  }
  
  chainSummary$alpha = chainSummary$alpha/ns
  
  chainSummary$meanvec = array(0, c(J, p, K))
  chainSummary$meanvec0 = matrix(0, p, K)
  for (k in 1:K) {
    del.om = transform_params(chainSummary$Sigma[,,k], chainSummary$alpha[,k])
    for (j in 1:J) {
      chainSummary$meanvec[j,,k] = chainSummary$xi[j,,k] + del.om$omega %*% del.om$delta*sqrt(2/pi)
    }
    chainSummary$meanvec0[,k] = chainSummary$xi0[,k] + del.om$omega %*% del.om$delta*sqrt(2/pi)
  }
  
  chainSummary$t = apply(chain$t,2,Mode)
  chainSummary$eta = mean(chain$eta)
  
  chainSummary
}