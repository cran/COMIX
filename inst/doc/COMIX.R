## ------------------------------------------------------------------------
knitr::opts_chunk$set(
  comment = '', fig.width = 6, fig.height = 6
)

## ------------------------------------------------------------------------
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
pmc <- list(naprt = 5, nburn = 200, nsave = 200)
# Fit the model:
res <- comix(Y, C, pmc = pmc, prior = prior)

# Relabel to resolve potential label switching issues:
res_relab <- relabelChain(res)

# Generate calibrated data:
cal <- calibrateNoDist(res_relab)

# Compare raw and calibrated data:
oldpar <- par(mfrow = c(1, 2)) 
plot(Y, col = C, xlim = range(Y[,1]), ylim = range(Y[,2]), main = "Raw Data",
     xlab = expression(Y[1]), ylab = expression(Y[2]))
plot(cal$Y_cal, col = C, main = "Calibrated Data",
     xlab = expression(Y[1]), ylab = expression(Y[2]))
par(oldpar)

# Get posterior estimates for the model parameters:
res_summary <- summarizeChain(res_relab)
# Check for instance, the cluster assignment labels:
table(res_summary$t)
# Indeed the same as 
colSums(njk)

# Or examine the skewness parameter for the non-trivial clusters:
res_summary$alpha[ , unique(res_summary$t)]
# And compare those to
cbind(alpha1, alpha2)

