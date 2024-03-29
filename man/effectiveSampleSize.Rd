% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effective_sample_size.R
\name{effectiveSampleSize}
\alias{effectiveSampleSize}
\title{This function creates an object that summarizes the effective sample size
for the parameters of the model.}
\usage{
effectiveSampleSize(res, params = c("w", "xi", "xi0", "psi", "G", "E", "eta"))
}
\arguments{
\item{res}{An object of class \code{COMIX} or \code{tidyChainCOMIX}.}

\item{params}{A character vector naming the parameters to compute the 
effective sample size for.}
}
\value{
An \code{effectiveSampleSizeCOMIX} object which is a named list,
with a named element for each requested parameter. Each element is a data
frame that includes the effective sample size for the parameter.
}
\description{
This function creates an object that summarizes the effective sample size
for the parameters of the model.
}
\examples{
library(COMIX)
# Number of observations for each sample (row) and cluster (column):
njk <- 
  matrix(
    c(
      150, 300,
      250, 200
    ),
    nrow = 2,
    byrow = TRUE
  )

# Dimension of data:
p <- 3

# Scale and skew parameters for first cluster:
Sigma1 <- matrix(0.5, nrow = p, ncol = p) + diag(0.5, nrow = p)
alpha1 <- rep(0, p)
alpha1[1] <- -5
# location parameter for first cluster in first sample:
xi11 <- rep(0, p)
# location parameter for first cluster in second sample (aligned with first):
xi21 <- rep(0, p)

# Scale and skew parameters for second cluster:
Sigma2 <- matrix(-1/3, nrow = p, ncol = p) + diag(1 + 1/3, nrow = p)
alpha2 <- rep(0, p)
alpha2[2] <- 5
# location parameter for second cluster in first sample:
xi12 <- rep(3, p)
# location parameter for second cluster in second sample (misaligned with first):
xi22 <- rep(4, p)

# Sample data:
set.seed(1)
Y <- 
  rbind(
    sn::rmsn(njk[1, 1], xi = xi11, Omega = Sigma1, alpha = alpha1),
    sn::rmsn(njk[1, 2], xi = xi12, Omega = Sigma2, alpha = alpha2),
    sn::rmsn(njk[2, 1], xi = xi21, Omega = Sigma1, alpha = alpha1),
    sn::rmsn(njk[2, 2], xi = xi22, Omega = Sigma2, alpha = alpha2)
  )

C <- c(rep(1, rowSums(njk)[1]), rep(2, rowSums(njk)[2]))

prior <- list(zeta = 1, K = 10)
pmc <- list(naprt = 5, nburn = 200, nsave = 200) # Reasonable usage
pmc <- list(naprt = 5, nburn = 2, nsave = 5) # Minimal usage for documentation
# Fit the model:
res <- comix(Y, C, pmc = pmc, prior = prior)

# Relabel to resolve potential label switching issues:
res_relab <- relabelChain(res)
effssz <- effectiveSampleSize(res_relab, "w")
# Or:
tidy_chain <- tidyChain(res_relab, "w")
effssz <- effectiveSampleSize(tidy_chain, "w")
# (see vignette for a more detailed example)
}
