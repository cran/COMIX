#' The function computes (and by default plots) estimates of the autocovariance
#' or autocorrelation function for the different parameters of the model. This
#' is a wrapper for coda::acf.
#' @param res An object of class \code{COMIX} or \code{tidyChainCOMIX}.
#' @param params A character vector naming the parameters to compute and plot 
#' the autocorrelation plots for.
#' @param only_non_trivial_clusters Logical, if \code{TRUE} only compute and/or
#' plot the autocorrelation for the clusters that are estimated to be non-empty.
#' @param lag.max maximum lag at which to calculate the autocorrelation. See more
#' details at ?acf.
#' @param type Character string giving the type of autocorrelation to be 
#' computed. See more details at ?acf.
#' @param plot Logical. If \code{TRUE} (the default) the autocorrelation is 
#' plotted.
#' @param ... Other arguments passed to \code{acf}.
#' @return An \code{acfParamsCOMIX} object which is a named list,
#' with a named element for each requested parameter. Each element is 
#' an object of class \code{acf} (from the \code{coda} package).
#' #' @examples
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
#' effssz <- effectiveSampleSize(res_relab, "w")
#' # Or:
#' tidy_chain <- tidyChain(res_relab, "w")
#' acf_w <- acfParams(tidy_chain, "w")
#' 
#' # (see vignette for a more detailed example)
#' @export
acfParams <- 
  function(
    res, 
    params = c("w", "xi", "xi0", "psi", "G", "E", "eta"), 
    only_non_trivial_clusters = TRUE,
    # coda::acf() parameters:
    lag.max = NULL, type = c("correlation", "covariance", "partial"),
    plot = TRUE, ...
    ) {
  stopifnot(is(res, "COMIX") | is(res, "tidyChainCOMIX"))
  
  if (is(res, "COMIX")) {
    tidy_chain <- tidyChain(res, params)
  } else {
    tidy_chain <- res
  }
  
  n <- attributes(tidy_chain)$n
  P <- attributes(tidy_chain)$p
  nsave <- attributes(tidy_chain)$nsave
  K <- attributes(tidy_chain)$K
  J <- attributes(tidy_chain)$J
  
  non_trivial_k <- attributes(tidy_chain)$non_trivial_k
  non_triv_j_k <- attributes(tidy_chain)$non_triv_j_k
  
  acfParams <- list()
  class(acfParams) <- "acfParamsCOMIX"
  attributes(acfParams)$n <- n
  attributes(acfParams)$p <- P
  attributes(acfParams)$nsave <- nsave
  attributes(acfParams)$K <- K
  attributes(acfParams)$J <- J
  attributes(acfParams)$non_trivial_k <- non_trivial_k
  attributes(acfParams)$non_triv_j_k <- non_triv_j_k

  # w -----
  if ("w" %in% params) {
    if (only_non_trivial_clusters) {
      a <- tidy_chain$w
      a$triv <- TRUE
      for (j in 1:J) {
        a$triv[a$j == j & a$k %in% non_triv_j_k[[as.character(j)]]] <- FALSE
      }
      a <- a %>% filter(!.data$triv) %>% select(-.data$triv)
    } else {
      a <- tidy_chain$w
    }
    
    a <- 
      a %>% 
      mutate(k = paste0("k=", .data$k), j = paste0("j=", .data$j)) %>%
      pivot_wider(names_from = c(.data$k, .data$j), values_from = .data$W) %>%
      select(-.data$iter)
    
    
    aa <- mcmc(data = a, start = 1)
    acfParams$w <-
      acf(x = aa, lag.max = lag.max, type = type, plot = plot, ...)
  }
  

  # xi0 -----
  if ("xi0" %in% params) {
    if (only_non_trivial_clusters) {
      a <- tidy_chain$xi0 %>% filter(.data$k %in% non_trivial_k)
    } else {
      a <- tidy_chain$xi0
    }
    
    a <-
      a %>%
      mutate(k = paste0("k=", .data$k), p = paste0("p=", .data$p)) %>%
      pivot_wider(names_from = c(.data$k, .data$p), values_from = c(.data$xi0)) %>%
      select(-.data$iter)

    aa <- mcmc(data = a, start = 1)
    acfParams$xi0 <- 
      acf(x = aa, lag.max = lag.max, type = type, plot = plot, ...)
  }
  

  # xi -----
  if ("xi" %in% params) {
    if (only_non_trivial_clusters) {
      a <- tidy_chain$xi
      a$triv <- TRUE
      for (j in 1:J) {
        a$triv[a$j == j & a$k %in% non_triv_j_k[[as.character(j)]]] <- FALSE
      }
      a <- a %>% filter(!.data$triv) %>% select(-.data$triv)
    } else {
      a <- tidy_chain$xi
    }
    
    a <-
      a %>%
      mutate(k = paste0("k=", .data$k), p = paste0("p=", .data$p), j = paste0("j=", .data$j)) %>%
      pivot_wider(names_from = c(.data$k, .data$p, .data$j), values_from = c(.data$xi)) %>%
      select(-.data$iter)

    aa <- mcmc(data = a, start = 1)
    acfParams$xi0 <- 
      acf(x = aa, lag.max = lag.max, type = type, plot = plot, ...)
  }

  # psi -----
  if ("psi" %in% params) {
    if (only_non_trivial_clusters) {
      a <- tidy_chain$psi %>% filter(.data$k %in% non_trivial_k)
    } else {
      a <- tidy_chain$psi
    }
    
    a <-
      a %>%
      mutate(k = paste0("k=", .data$k), p = paste0("p=", .data$p)) %>%
      pivot_wider(names_from = c(.data$k, .data$p), values_from = c(.data$psi)) %>%
      select(-.data$iter)

    aa <- mcmc(data = a, start = 1)
    acfParams$psi <- 
      acf(x = aa, lag.max = lag.max, type = type, plot = plot, ...)
  }

  # G -----
  if ("G" %in% params) {
    if (only_non_trivial_clusters) {
      a <- tidy_chain$G %>% filter(.data$k %in% non_trivial_k)
    } else {
      a <- tidy_chain$G
    }
    
    a <-
      a %>%
      mutate(k = paste0("k=", .data$k), p1 = paste0("p1=", .data$p1), p2 = paste0("p2=", .data$p2)) %>%
      pivot_wider(names_from = c(.data$k, .data$p1, .data$p2), values_from = c(.data$G)) %>%
      select(-.data$iter)

    aa <- mcmc(data = a, start = 1)
    acfParams$G <- 
      acf(x = aa, lag.max = lag.max, type = type, plot = plot, ...)  
    }

  # E -----
  if ("E" %in% params) {
    if (only_non_trivial_clusters) {
      a <- tidy_chain$E %>% filter(.data$k %in% non_trivial_k)
    } else {
      a <- tidy_chain$E
    }
    
    a <-
      a %>%
      mutate(k = paste0("k=", .data$k), p1 = paste0("p1=", .data$p1), p2 = paste0("p2=", .data$p2)) %>%
      pivot_wider(names_from = c(.data$k, .data$p1, .data$p2), values_from = c(.data$E)) %>%
      select(-.data$iter)
    
    aa <- mcmc(data = a, start = 1)
    acfParams$E <- 
      acf(x = aa, lag.max = lag.max, type = type, plot = plot, ...)  
  }


  # eta -----
  if ("eta" %in% params) {
    eta <- mcmc(data = pull(tidy_chain$eta), start = 1)
    acfParams$eta <- 
      acf(x = eta, lag.max = lag.max, type = type, plot = plot, ...)
  }
  
  return(acfParams)
}
