#' This function creates trace plots for different parameters of the MCMC chain.
#'
#' @param res An object of class \code{COMIX} or \code{tidyChainCOMIX}.
#' @param param Character, naming the parameter to create a trace plot for.
#' @return A \code{ggplot2} plot containing the trace plot.
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
#' plotTracePlots(res_relab, "w")
#' # Or:
#' tidy_chain <- tidyChain(res_relab, "w")
#' plotTracePlots(tidy_chain, "w")
#' # (see vignette for a more detailed example)
#' @export
plotTracePlots <- function(res, param) {
  stopifnot(is(res, "COMIX") | is(res, "tidyChainCOMIX"))
  stopifnot(length(param) == 1)
  stopifnot(param %in% c("w", "xi0", "xi", "psi", "G", "E", "eta"))
  
  if (is(res, "COMIX")) {
    tidy_chain <- tidyChain(res, param)
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
  
  glob_freq_t <- attributes(tidy_chain)$glob_freq_t

  frq_t <- glob_freq_t$frq_t[glob_freq_t$k %in% non_trivial_k]
  k_names <- paste0("Cluster ", non_trivial_k, "\n(Est. Freq. = ", round(frq_t, 2), ")")
  names(k_names) <- non_trivial_k

  j_names <- paste0("Sample ", 1:J)
  names(j_names) <- 1:J
  
  
  # w -----
  if (param == "w") {
    g <- 
      tidy_chain$w %>% 
      ggplot() +
      geom_point(
        mapping = aes(x = .data$iter, y = .data$W, color = .data$k), 
        show.legend = FALSE, 
        size = 1,
        alpha = 0.4) +
      facet_wrap(~ .data$j, labeller = labeller(j = j_names)) +
      ylab(expression(omega)) +
      xlab("Iteration") +
      theme(axis.text.x = element_text(angle = -90, vjust = 0.5))

    return(g)
  }
  
  # xi0 -----
  if (param == "xi0") {
    g <- 
      tidy_chain$xi0 %>%
      filter(.data$k %in% non_trivial_k) %>%
      ggplot() +
      geom_point(
        mapping = aes(x = .data$iter, y = .data$xi0, color = .data$p), 
        show.legend = FALSE, 
        size = 1,
        alpha = 0.4
        ) +
      facet_wrap( ~ .data$k, labeller = labeller(k = k_names)) +
      ylab(expression(xi[0])) +
      xlab("Iteration") +
      theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
  
    return(g)
  }
  
  
  # xi -----
  if (param == "xi") {
    a <- tidy_chain$xi
    a$triv <- TRUE

    for (j in 1:J) {
      a$triv[a$j == j & a$k %in% non_triv_j_k[[as.character(j)]]] <- FALSE
    }
    
    g <-
      a %>% 
      filter(!.data$triv) %>%
      ggplot() +
      geom_point(
        mapping = aes(x = .data$iter, y = .data$xi, color = .data$p), 
        show.legend = FALSE, 
        size = 1,
        alpha = 0.4
        ) +
      facet_grid(.data$j ~ .data$k, labeller = labeller(j = j_names, k = k_names)) +
      ylab(expression(xi)) +
      xlab("Iteration") +
      theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
    
    return(g)
  }
  

  # psi -----
  if (param == "psi") {
    g <-
      tidy_chain$psi %>% 
      filter(.data$k %in% non_trivial_k) %>%
      ggplot() +
      geom_point(
        mapping = aes(x = .data$iter, y = .data$psi, color = .data$p), 
        show.legend = FALSE, 
        size = 1,
        alpha = 0.4
        ) +
      facet_wrap(~ .data$k, labeller = labeller(k = k_names)) +
      ylab(expression(psi)) +
      xlab("Iteration") +
      theme(axis.text.x = element_text(angle = -90, vjust = 0.5))
    
    return(g)
  }

    
  # G -----
  if ("G" %in% param) {
    g <-  
      tidy_chain$G %>%
      filter(.data$k %in% non_trivial_k) %>%
      mutate(pp = paste0(.data$p1, ",", .data$p2)) %>%
      ggplot() +
      geom_point(
        mapping = aes(x = .data$iter, y = .data$G, color = .data$pp), 
        show.legend = FALSE, 
        size = 1,
        alpha = 0.4
        ) +
      facet_wrap(~ .data$k, labeller = labeller(k = k_names)) +
      ylab("G") +
      xlab("Iteration") +
      theme(axis.text.x = element_text(angle = -90, vjust = 0.5))

    return(g)
  }
  
  # E -----
  if ("E" %in% param) {
    g <-    
      tidy_chain$E %>% 
      filter(.data$k %in% non_trivial_k) %>%
      mutate(pp = paste0(.data$p1, ",", .data$p2)) %>%
      ggplot() +
      geom_point(
        mapping = aes(x = .data$iter, y = .data$E, color = .data$pp), 
        show.legend = FALSE, 
        size = 1,
        alpha = 0.4
        ) +
      facet_wrap(~ .data$k, labeller = labeller(k = k_names)) +
      ylab("E") +
      xlab("Iteration") +
      theme(axis.text.x = element_text(angle = -90, vjust = 0.5))

    return(g)
  }
    
  # eta -----
  if ("eta" %in% param) {
    g <- 
      ggplot(tidy_chain$eta) + 
        geom_point(aes(x = .data$iter, y = .data$eta), size = 1, alpha = 0.4) +
        ylab(expression(eta)) +
        xlab("Iteration")
    return(g)
  }
  
  return(NULL)
}
