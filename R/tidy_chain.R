#' This function creates tidy versions of the stored chain. This object can then be used as input 
#' for the other diagnostic functions in this package.
#'
#' @param res An object of class COMIX.
#' @param params A character vector naming the parameters to tidy.
#' @return A \code{tidyChainCOMIX} object: a named list of class whose length is the length 
#' of \code{params}. Each element of the list contains a tibble with a tidy version of the samples
#' from the MCMC chain.
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
#' tidy_chain <- tidyChain(res_relab)
#' # (see vignette for a more detailed example)
#' @export
tidyChain <- function(res, params = c("t", "w", "xi", "xi0", "psi", "G", "E", "eta", "Sigma", "alpha")) {
  stopifnot(is(res, "COMIX"))
  
  sum_chain <- COMIX::summarizeChain(res)
  
  n <- nrow(res$data$Y)
  P <- ncol(res$data$Y)
  nsave <- res$pmc$nsave
  K <- res$prior$K
  J <- length(unique(res$data$C))
  
  non_trivial_k <- sort(unique(sum_chain$t))
  non_triv_j_k <- split(x = sum_chain$t, f = res$data$C) %>% lapply(function(t) sort(unique(t)))
  
  tidy_chain <- list()
  class(tidy_chain) <- "tidyChainCOMIX"
  attributes(tidy_chain)$pmc <- res$pmc
  attributes(tidy_chain)$prior <- res$prior
  attributes(tidy_chain)$n <- n
  attributes(tidy_chain)$p <- P
  attributes(tidy_chain)$nsave <- nsave
  attributes(tidy_chain)$K <- K
  attributes(tidy_chain)$J <- J
  attributes(tidy_chain)$non_trivial_k <- non_trivial_k
  attributes(tidy_chain)$non_triv_j_k <- non_triv_j_k
  attributes(tidy_chain)$glob_freq_t <- 
    data.frame(
      x = (P + 1) / 2, 
      table(sum_chain$t) / sum(table(sum_chain$t))
    ) %>%
    rename(c("k" = "Var1", "frq_t" = "Freq"))
  attributes(tidy_chain)$local_freq_t <-
    tibble(x = (P + 1) / 2, k = sum_chain$t, j = c(res$data$C)) %>%
    group_by(.data$j) %>%
    # summarize(Freq = table(.data$k) / sum(table(.data$k))) %>%
    reframe(Freq = table(.data$k) / sum(table(.data$k))) %>%
    mutate(x = (P + 1) / 2, k = names(.data$Freq), frq_t = as.numeric(.data$Freq)) %>%
    select(-.data$Freq)


  # t -----
  if ("t" %in% params) {
    tidy_chain$t <- tibble(t = res$chain$t)
  }
  
  # w -----
  if ("w" %in% params) {
    w_tb <-
      expand_grid(j = 1:J, k = 1:K) %>%
      slice(rep(1:n(), each = nsave)) %>%
      mutate(
        j = factor(.data$j),
        k = factor(.data$k),
        W = 0
      )
    w_tb$iter <- rep(1:nsave, J * K)
    
    for ( j in 1:J ) {
      for ( k in 1:K ) {
        w_tb$W[w_tb$j == j & w_tb$k == k] <- res$chain$W[j, k, ]
      }
    }
    
    tidy_chain$w <- w_tb
  }
  
  
  # xi0 -----
  if ("xi0" %in% params) {
    xi0_tb <-
      expand_grid(p = 1:P, k = 1:K) %>%
      slice(rep(1:n(), each = nsave)) %>%
      mutate(
        p = factor(.data$p),
        k = factor(.data$k),
        xi0 = 0
      )
    xi0_tb$iter <- rep(1:nsave, P * K)
    
    for ( p in 1:P ) {
      for ( k in 1:K ) {
        xi0_tb$xi0[xi0_tb$p == p & xi0_tb$k == k] <- res$chain$xi0[p, k, ]
      }
    }
    
    tidy_chain$xi0 <- xi0_tb
  }
  
  # xi -----
  if ("xi" %in% params) {
    xi_tb <- 
      tibble(
        k = NA_integer_, 
        p = NA_integer_, 
        j = NA_integer_, 
        xi = NA_real_, 
        iter = rep(1:nsave, K * P * J)
      )
    
    counter <- 1
    for ( k in 1:K ) {
      for ( p in 1:P ) {
        for ( j in 1:J ) {
          current_idx <- counter : (counter + nsave - 1)
          xi_tb$k[current_idx] <- k
          xi_tb$p[current_idx] <- p
          xi_tb$j[current_idx] <- j
          xi_tb$xi[current_idx] <- res$chain$xi[j, P * (k - 1) + p, ]
          counter <- counter + nsave
        }
      }
    }
    xi_tb <- mutate(xi_tb, p = factor(p), k = factor(k), j = factor(j))

    tidy_chain$xi <- xi_tb    
  }
  
  # psi -----
  if ("psi" %in% params) {
    psi_tb <-
      expand_grid(p = 1:P, k = 1:K) %>%
      slice(rep(1:n(), each = nsave)) %>%
      mutate(
        p = factor(.data$p),
        k = factor(.data$k),
        psi = 0
      )
    psi_tb$iter <- rep(1:nsave, P * K)
    
    for ( p in 1:P ) {
      for ( k in 1:K ) {
        psi_tb$psi[psi_tb$p == p & psi_tb$k == k] <- res$chain$psi[p, k, ]
      }
    }
   
    tidy_chain$psi <- psi_tb
  }
  
  # G -----
  if ("G" %in% params) {
    G_tb <- 
      tibble(
        k = NA_integer_, 
        p1 = NA_integer_, 
        p2 = NA_integer_, 
        G = NA_real_, 
        iter = rep(1:nsave, K * (P + 1) * P / 2)
      )
    
    counter <- 1
    for ( k in 1:K ) {
      for ( p1 in 1:P ) {
        for ( p2 in p1:P ) {
          current_idx <- counter : (counter + nsave - 1)
          G_tb$k[current_idx] <- k
          G_tb$p1[current_idx] <- p1
          G_tb$p2[current_idx] <- p2
          G_tb$G[current_idx] <- res$chain$G[p1, P * (k - 1) + p2, ]
          counter <- counter + nsave
        }
      }
    }
    G_tb <- mutate(G_tb, k = factor(k), p1 = factor(p1), p2 = factor(p2))

    tidy_chain$G <- G_tb
  }
  
  # E -----
  if ("E" %in% params) {
    E_tb <- 
      tibble(
        k = NA_integer_, 
        p1 = NA_integer_, 
        p2 = NA_integer_, 
        E = NA_real_, 
        iter = rep(1:nsave, K * (P + 1) * P / 2)
      )
    
    counter <- 1
    for ( k in 1:K ) {
      for ( p1 in 1:P ) {
        for ( p2 in p1:P ) {
          current_idx <- counter : (counter + nsave - 1)
          E_tb$k[current_idx] <- k
          E_tb$p1[current_idx] <- p1
          E_tb$p2[current_idx] <- p2
          E_tb$E[current_idx] <- res$chain$E[p1, P * (k - 1) + p2, ]
          counter <- counter + nsave
        }
      }
    }
    E_tb <- mutate(E_tb, k = factor(k), p1 = factor(p1), p2 = factor(p2))

    tidy_chain$E <- E_tb
  }
  
  # eta -----
  if ("eta" %in% params) {
    tidy_chain$eta <- tibble(eta = res$chain$eta) %>% mutate(iter = row_number())
  }

  # Sigma -----
  if ("Sigma" %in% params) {
    Sigma_tb <- 
      tibble(
        k = NA_integer_, 
        p1 = NA_integer_, 
        p2 = NA_integer_, 
        Sigma = NA_real_, 
        iter = rep(1:nsave, K * (P + 1) * P / 2)
      )
    
    counter <- 1
    for ( k in 1:K ) {
      for ( p1 in 1:P ) {
        for ( p2 in p1:P ) {
          current_idx <- counter : (counter + nsave - 1)
          Sigma_tb$k[current_idx] <- k
          Sigma_tb$p1[current_idx] <- p1
          Sigma_tb$p2[current_idx] <- p2
          Sigma_tb$Sigma[current_idx] <- res$chain$Sigma[p1, P * (k - 1) + p2, ]
          counter <- counter + nsave
        }
      }
    }
    Sigma_tb <- mutate(Sigma_tb, k = factor(k), p1 = factor(p1), p2 = factor(p2))
    
    tidy_chain$Sigma <- Sigma_tb
  }
  
  # alpha -----
  if ("alpha" %in% params) {
    alpha_tb <-
      expand_grid(p = 1:P, k = 1:K) %>%
      slice(rep(1:n(), each = nsave)) %>%
      mutate(
        p = factor(.data$p),
        k = factor(.data$k),
        alpha = 0
      )
    alpha_tb$iter <- rep(1:nsave, P * K)
    
    for ( p in 1:P ) {
      for ( k in 1:K ) {
        alpha_tb$alpha[alpha_tb$p == p & alpha_tb$k == k] <- res$chain$alpha[p, k, ]
      }
    }
    
    tidy_chain$alpha <- alpha_tb
  }
  
  return(tidy_chain)
}
